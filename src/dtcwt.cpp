#include "dtcwt.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <string>

namespace fs = std::filesystem;

namespace {

inline int reflect_index(int i, int n) {
    // Symmetric extension (reflection), consistent for forward/inverse
    // ... 2 1 0 1 2 3 2 1 ...
    if (n <= 1) return 0;
    while (i < 0 || i >= n) {
        if (i < 0) i = -i - 1;
        if (i >= n) i = 2 * n - i - 1;
    }
    return i;
}

static void check_real_image(const RealImage2D& im) {
    if (im.width <= 0 || im.height <= 0) {
        throw std::runtime_error("RealImage2D: invalid dimensions");
    }
    if (im.data.size() != static_cast<size_t>(im.width) * static_cast<size_t>(im.height)) {
        throw std::runtime_error("RealImage2D: data size mismatch");
    }
}

static RealImage2D make_real(int w, int h) {
    RealImage2D out{w, h, {}};
    out.data.assign(static_cast<size_t>(w) * static_cast<size_t>(h), 0.0);
    return out;
}

static ComplexImage2D make_cplx(int w, int h) {
    ComplexImage2D out{w, h, {}};
    out.data.assign(static_cast<size_t>(w) * static_cast<size_t>(h), {0.0, 0.0});
    return out;
}

inline double& at(RealImage2D& im, int x, int y) {
    return im.data[static_cast<size_t>(y) * static_cast<size_t>(im.width) + static_cast<size_t>(x)];
}
inline double at(const RealImage2D& im, int x, int y) {
    return im.data[static_cast<size_t>(y) * static_cast<size_t>(im.width) + static_cast<size_t>(x)];
}

inline std::complex<double>& at(ComplexImage2D& im, int x, int y) {
    return im.data[static_cast<size_t>(y) * static_cast<size_t>(im.width) + static_cast<size_t>(x)];
}

// -------------------------------
// Filter definitions (STRICT per AI_DTCWT_Qshift_Spec.md)
// - Level 1: near-symmetric (biort) 5/7 (near_sym_a)
// - Level >= 2: Kingsbury Q-shift (qshift_a) length 10
//
// NOTE:
// - Tree A uses shift=0
// - Tree B uses shift=1
// - Tree B filters are time-reversed counterparts (for qshift provided explicitly;
//   for near-symmetric we derive Tree B by reversing the provided arrays).
// -------------------------------

struct FilterBank1D {
    // Analysis
    std::vector<double> h0a, h1a, h0b, h1b;
    // Synthesis
    std::vector<double> g0a, g1a, g0b, g1b;
};

static std::vector<double> reversed(const std::vector<double>& v) {
    std::vector<double> r = v;
    std::reverse(r.begin(), r.end());
    return r;
}

// Standard biorthogonal QMF: reverse + alternating sign
static std::vector<double> qmf(const std::vector<double>& f) {
    std::vector<double> r(f.rbegin(), f.rend());
    for (size_t i = 0; i < r.size(); ++i) {
        if (i % 2 == 1) r[i] = -r[i];
    }
    return r;
}

static FilterBank1D near_sym_a_filterbank() {
    FilterBank1D fb;
    // From: coeffs.biort('near_sym_a')
    // h0a: Tree A analysis lowpass
    // g0a: Tree A synthesis lowpass
    // h1a: Tree A analysis highpass
    // g1a: Tree A synthesis highpass
    // h0b: Tree B analysis lowpass
    // g0b: Tree B synthesis lowpass
    // h1b: Tree B analysis highpass
    // g1b: Tree B synthesis highpass
    
    // Tree A
    fb.h0a = { -0.050000000000000003, 0.25, 0.59999999999999998, 0.25, -0.050000000000000003 };
    fb.g0a = { -0.010714285714285713, -0.053571428571428568, 0.26071428571428573, 0.6071428571428571, 0.26071428571428573, -0.053571428571428568, -0.010714285714285713 };
    fb.h1a = { -0.010714285714285713, 0.053571428571428568, 0.26071428571428573, -0.6071428571428571, 0.26071428571428573, 0.053571428571428568, -0.010714285714285713 };
    fb.g1a = { -0.050000000000000003, -0.25, 0.59999999999999998, -0.25, -0.050000000000000003 };

    // Tree B
    fb.h0b = { 0.010714285714285713, -0.053571428571428568, -0.26071428571428573, 0.6071428571428571, -0.26071428571428573, -0.053571428571428568, 0.010714285714285713 };
    fb.g0b = { -0.050000000000000003, -0.25, 0.59999999999999998, -0.25, -0.050000000000000003 };
    fb.h1b = { -0.050000000000000003, 0.25, 0.59999999999999998, 0.25, -0.050000000000000003 };
    fb.g1b = { 0.010714285714285713, 0.053571428571428568, -0.26071428571428573, -0.6071428571428571, -0.26071428571428573, 0.053571428571428568, 0.010714285714285713 };

    return fb;
}

static FilterBank1D qshift10_filterbank() {
    FilterBank1D fb;
    // From: coeffs.qshift('qshift_a')
    // Tree A
    fb.h0a = { 0.051130405283831656, -0.013975370246888838, -0.10983605166597087, 0.26383956105893763, 0.76662846779303717, 0.56365571012705151, 0.00087362269521709679, -0.1002312195074762, -0.0016896812725281543, -0.0061818818921164382 };
    fb.g0a = { -0.0061818818921164382, -0.0016896812725281543, -0.1002312195074762, 0.00087362269521709679, 0.56365571012705151, 0.76662846779303717, 0.26383956105893763, -0.10983605166597087, -0.013975370246888838, 0.051130405283831656 };
    fb.h1a = { -0.0061818818921164382, 0.0016896812725281543, -0.1002312195074762, -0.00087362269521709679, 0.56365571012705151, -0.76662846779303717, 0.26383956105893763, 0.10983605166597087, -0.013975370246888838, -0.051130405283831656 };
    fb.g1a = { -0.051130405283831656, -0.013975370246888838, 0.10983605166597087, 0.26383956105893763, -0.76662846779303717, 0.56365571012705151, -0.00087362269521709679, -0.1002312195074762, 0.0016896812725281543, -0.0061818818921164382 };

    // Tree B
    fb.h0b = { -0.0061818818921164382, -0.0016896812725281543, -0.1002312195074762, 0.00087362269521709679, 0.56365571012705151, 0.76662846779303717, 0.26383956105893763, -0.10983605166597087, -0.013975370246888838, 0.051130405283831656 };
    fb.g0b = { 0.051130405283831656, -0.013975370246888838, -0.10983605166597087, 0.26383956105893763, 0.76662846779303717, 0.56365571012705151, 0.00087362269521709679, -0.1002312195074762, -0.0016896812725281543, -0.0061818818921164382 };
    fb.h1b = { -0.051130405283831656, -0.013975370246888838, 0.10983605166597087, 0.26383956105893763, -0.76662846779303717, 0.56365571012705151, -0.00087362269521709679, -0.1002312195074762, 0.0016896812725281543, -0.0061818818921164382 };
    fb.g1b = { -0.0061818818921164382, 0.0016896812725281543, -0.1002312195074762, -0.00087362269521709679, 0.56365571012705151, -0.76662846779303717, 0.26383956105893763, 0.10983605166597087, -0.013975370246888838, -0.051130405283831656 };
    return fb;
}

// -------------------------------
// 1D analysis/synthesis (reflection boundary)
// shift=0 for Tree A, shift=1 for Tree B (per spec)
// -------------------------------

static void analysis_1d(const double* x, int n,
                        const std::vector<double>& h0,
                        const std::vector<double>& h1,
                        int shift,
                        std::vector<double>& lo,
                        std::vector<double>& hi)
{
    const int L0 = static_cast<int>(h0.size());
    const int L1 = static_cast<int>(h1.size());
    if (L0 <= 0 || L1 <= 0) {
        throw std::runtime_error("analysis_1d: invalid filter lengths");
    }
    const int m = (n + 1) / 2; // ceil(n/2)
    lo.assign(m, 0.0);
    hi.assign(m, 0.0);

    for (int k = 0; k < m; ++k) {
        const int center = 2 * k + shift;
        double acc0 = 0.0;
        double acc1 = 0.0;
        for (int t = 0; t < L0; ++t) {
            const int idx = reflect_index(center - t, n);
            acc0 += h0[static_cast<size_t>(t)] * x[idx];
        }
        for (int t = 0; t < L1; ++t) {
            const int idx = reflect_index(center - t, n);
            acc1 += h1[static_cast<size_t>(t)] * x[idx];
        }
        lo[k] = acc0;
        hi[k] = acc1;
    }
}

// Q-shift analysis 1D (coldfilt)
// x        : input signal, length n
// h0, h1   : even-length Q-shift analysis filters (e.g. h0a/h1a or h0b/h1b)
// lo, hi   : outputs, length = n/2
// symmetric extension is used
static void analysis_1d_qshift_coldfilt(
    const std::vector<double>& x,
    const std::vector<double>& h0,
    const std::vector<double>& h1,
    int shift,
    std::vector<double>& lo,
    std::vector<double>& hi)
{
    // NOTE:
    // We treat Q-shift stages as a critically-sampled 2-channel filterbank with
    // explicit phase control (shift = 0/1). This must match synthesis exactly.
    const int n = static_cast<int>(x.size());
    const int L0 = static_cast<int>(h0.size());
    const int L1 = static_cast<int>(h1.size());
    if ((n % 2) != 0) throw std::runtime_error("analysis_1d_qshift_coldfilt: input length must be even");
    if (L0 <= 0 || L1 <= 0) throw std::runtime_error("analysis_1d_qshift_coldfilt: invalid filter lengths");

    const int m = n / 2;
    lo.assign(m, 0.0);
    hi.assign(m, 0.0);

    // Phase-aligned filtering: center = 2*k + shift
    for (int k = 0; k < m; ++k) {
        const int center = 2 * k + shift;
        double acc0 = 0.0;
        double acc1 = 0.0;

        for (int t = 0; t < L0; ++t) {
            const int idx = reflect_index(center - t, n);
            acc0 += h0[static_cast<size_t>(t)] * x[static_cast<size_t>(idx)];
        }
        for (int t = 0; t < L1; ++t) {
            const int idx = reflect_index(center - t, n);
            acc1 += h1[static_cast<size_t>(t)] * x[static_cast<size_t>(idx)];
        }

        lo[static_cast<size_t>(k)] = acc0;
        hi[static_cast<size_t>(k)] = acc1;
    }
}

static void synth_1d(const std::vector<double>& lo,
                     const std::vector<double>& hi,
                     int n,
                     const std::vector<double>& g0,
                     const std::vector<double>& g1,
                     int shift,
                     std::vector<double>& y)
{
    const int L0 = static_cast<int>(g0.size());
    const int L1 = static_cast<int>(g1.size());
    if (L0 <= 0 || L1 <= 0) {
        throw std::runtime_error("synth_1d: invalid filter lengths");
    }
    if (lo.size() != hi.size()) {
        throw std::runtime_error("synth_1d: lo/hi size mismatch");
    }

    y.assign(static_cast<size_t>(n), 0.0);
    const int m = static_cast<int>(lo.size());
    for (int k = 0; k < m; ++k) {
        const int center = 2 * k + shift;
        for (int t = 0; t < L0; ++t) {
            // Synthesis corresponds to conv(upsample(lo,2), g0):
            // contribution of lo[k] is placed at output index (2*k + shift + t)
            const int idx = reflect_index(center + t, n);
            y[static_cast<size_t>(idx)] += g0[static_cast<size_t>(t)] * lo[static_cast<size_t>(k)];
        }
        for (int t = 0; t < L1; ++t) {
            const int idx = reflect_index(center + t, n);
            y[static_cast<size_t>(idx)] += g1[static_cast<size_t>(t)] * hi[static_cast<size_t>(k)];
        }
    }
    // IMPORTANT: compensate for upsampling-by-2
    for (double& v : y) v *= 2.0;
}

// Q-shift synthesis 1D (colifilt)
// lo, hi  : low/high subbands, length m
// g0, g1  : even-length Q-shift synthesis filters
// x       : reconstructed signal, length = 2*m
static void synth_1d_qshift_colifilt(
    const std::vector<double>& lo,
    const std::vector<double>& hi,
    const std::vector<double>& g0,
    const std::vector<double>& g1,
    int shift,
    std::vector<double>& x)
{
    // NOTE:
    // Must be the exact dual of analysis_1d_qshift_coldfilt() with the same shift.
    if (lo.size() != hi.size()) throw std::runtime_error("synth_1d_qshift_colifilt: lo/hi size mismatch");
    const int m = static_cast<int>(lo.size());
    const int n = 2 * m;
    const int L0 = static_cast<int>(g0.size());
    const int L1 = static_cast<int>(g1.size());
    if ((n % 2) != 0) throw std::runtime_error("synth_1d_qshift_colifilt: output length must be even");
    if (L0 <= 0 || L1 <= 0) throw std::runtime_error("synth_1d_qshift_colifilt: invalid filter lengths");

    x.assign(static_cast<size_t>(n), 0.0);

    for (int k = 0; k < m; ++k) {
        const int center = 2 * k + shift;
        for (int t = 0; t < L0; ++t) {
            const int idx = reflect_index(center + t, n);
            x[static_cast<size_t>(idx)] += g0[static_cast<size_t>(t)] * lo[static_cast<size_t>(k)];
        }
        for (int t = 0; t < L1; ++t) {
            const int idx = reflect_index(center + t, n);
            x[static_cast<size_t>(idx)] += g1[static_cast<size_t>(t)] * hi[static_cast<size_t>(k)];
        }
    }
}

// -------------------------------
// 2D separable analysis for a (rowTree, colTree)
// Row analysis then column analysis (strict order per spec)
// -------------------------------

struct Subbands2D {
    RealImage2D LL, LH, HL, HH;
};

static Subbands2D analysis_2d_tree(const RealImage2D& in,
                                  const std::vector<double>& row_h0,
                                  const std::vector<double>& row_h1,
                                  int row_shift,
                                  const std::vector<double>& col_h0,
                                  const std::vector<double>& col_h1,
                                  int col_shift)
{
    const int W = in.width;
    const int H = in.height;
    const bool row_is_qshift = (row_h0.size() % 2 == 0) && (row_h1.size() % 2 == 0);
    const bool col_is_qshift = (col_h0.size() % 2 == 0) && (col_h1.size() % 2 == 0);

    if (row_is_qshift && (W % 2 != 0)) {
        throw std::runtime_error("analysis_2d_tree: Q-shift row filtering requires even width");
    }
    if (col_is_qshift && (H % 2 != 0)) {
        throw std::runtime_error("analysis_2d_tree: Q-shift col filtering requires even height");
    }

    const int W2 = row_is_qshift ? (W / 2) : ((W + 1) / 2);
    const int H2 = col_is_qshift ? (H / 2) : ((H + 1) / 2);

    // Row analysis
    RealImage2D Lr = make_real(W2, H);
    RealImage2D Hr = make_real(W2, H);
    std::vector<double> row(static_cast<size_t>(W));
    std::vector<double> lo, hi;

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) row[static_cast<size_t>(x)] = at(in, x, y);
        if (row_is_qshift) {
            analysis_1d_qshift_coldfilt(row, row_h0, row_h1, row_shift, lo, hi);
        } else {
            analysis_1d(row.data(), W, row_h0, row_h1, row_shift, lo, hi);
        }
        for (int x = 0; x < W2; ++x) {
            at(Lr, x, y) = lo[static_cast<size_t>(x)];
            at(Hr, x, y) = hi[static_cast<size_t>(x)];
        }
    }

    // Column analysis
    Subbands2D sb{make_real(W2, H2), make_real(W2, H2), make_real(W2, H2), make_real(W2, H2)};
    std::vector<double> colL(static_cast<size_t>(H));
    std::vector<double> colH(static_cast<size_t>(H));
    std::vector<double> lo_c, hi_c;

    for (int x = 0; x < W2; ++x) {
        for (int y = 0; y < H; ++y) {
            colL[static_cast<size_t>(y)] = at(Lr, x, y);
            colH[static_cast<size_t>(y)] = at(Hr, x, y);
        }
        if (col_is_qshift) {
            analysis_1d_qshift_coldfilt(colL, col_h0, col_h1, col_shift, lo_c, hi_c);
        } else {
            analysis_1d(colL.data(), H, col_h0, col_h1, col_shift, lo_c, hi_c);
        }
        for (int y = 0; y < H2; ++y) {
            at(sb.LL, x, y) = lo_c[static_cast<size_t>(y)];
            at(sb.LH, x, y) = hi_c[static_cast<size_t>(y)];
        }
        if (col_is_qshift) {
            analysis_1d_qshift_coldfilt(colH, col_h0, col_h1, col_shift, lo_c, hi_c);
        } else {
            analysis_1d(colH.data(), H, col_h0, col_h1, col_shift, lo_c, hi_c);
        }
        for (int y = 0; y < H2; ++y) {
            at(sb.HL, x, y) = lo_c[static_cast<size_t>(y)];
            at(sb.HH, x, y) = hi_c[static_cast<size_t>(y)];
        }
    }
    return sb;
}

static RealImage2D synth_2d_tree(const RealImage2D& LL,
                                 const RealImage2D& LH,
                                 const RealImage2D& HL,
                                 const RealImage2D& HH,
                                 int outW,
                                 int outH,
                                 const std::vector<double>& row_g0,
                                 const std::vector<double>& row_g1,
                                 int row_shift,
                                 const std::vector<double>& col_g0,
                                 const std::vector<double>& col_g1,
                                 int col_shift)
{
    // Column synthesis first (strict)
    const int W2 = LL.width;
    const int H2 = LL.height;
    if (LH.width != W2 || LH.height != H2 || HL.width != W2 || HL.height != H2 || HH.width != W2 || HH.height != H2) {
        throw std::runtime_error("synth_2d_tree: subband size mismatch");
    }

    const bool row_is_qshift = (row_g0.size() % 2 == 0) && (row_g1.size() % 2 == 0);
    const bool col_is_qshift = (col_g0.size() % 2 == 0) && (col_g1.size() % 2 == 0);
    if (col_is_qshift && (outH % 2 != 0)) {
        throw std::runtime_error("synth_2d_tree: Q-shift col synthesis requires even outH");
    }
    if (row_is_qshift && (outW % 2 != 0)) {
        throw std::runtime_error("synth_2d_tree: Q-shift row synthesis requires even outW");
    }

    RealImage2D Lr = make_real(W2, outH);
    RealImage2D Hr = make_real(W2, outH);

    std::vector<double> lo(H2), hi(H2), ycol;
    for (int x = 0; x < W2; ++x) {
        for (int y = 0; y < H2; ++y) {
            lo[static_cast<size_t>(y)] = at(LL, x, y);
            hi[static_cast<size_t>(y)] = at(LH, x, y);
        }
        if (col_is_qshift) {
            synth_1d_qshift_colifilt(lo, hi, col_g0, col_g1, col_shift, ycol);
        } else {
            synth_1d(lo, hi, outH, col_g0, col_g1, col_shift, ycol);
        }
        for (int y = 0; y < outH; ++y) at(Lr, x, y) = ycol[static_cast<size_t>(y)];

        for (int y = 0; y < H2; ++y) {
            lo[static_cast<size_t>(y)] = at(HL, x, y);
            hi[static_cast<size_t>(y)] = at(HH, x, y);
        }
        if (col_is_qshift) {
            synth_1d_qshift_colifilt(lo, hi, col_g0, col_g1, col_shift, ycol);
        } else {
            synth_1d(lo, hi, outH, col_g0, col_g1, col_shift, ycol);
        }
        for (int y = 0; y < outH; ++y) at(Hr, x, y) = ycol[static_cast<size_t>(y)];
    }

    // Row synthesis
    RealImage2D out = make_real(outW, outH);
    std::vector<double> rowL(static_cast<size_t>(W2));
    std::vector<double> rowH(static_cast<size_t>(W2));
    std::vector<double> yrow;
    for (int y = 0; y < outH; ++y) {
        for (int x = 0; x < W2; ++x) {
            rowL[static_cast<size_t>(x)] = at(Lr, x, y);
            rowH[static_cast<size_t>(x)] = at(Hr, x, y);
        }
        if (row_is_qshift) {
            synth_1d_qshift_colifilt(rowL, rowH, row_g0, row_g1, row_shift, yrow);
        } else {
            synth_1d(rowL, rowH, outW, row_g0, row_g1, row_shift, yrow);
        }
        for (int x = 0; x < outW; ++x) at(out, x, y) = yrow[static_cast<size_t>(x)];
    }

    return out;
}

static ComplexImage2D pack_complex(const RealImage2D& re, const RealImage2D& im) {
    if (re.width != im.width || re.height != im.height) throw std::runtime_error("pack_complex size mismatch");
    ComplexImage2D out = make_cplx(re.width, re.height);
    for (size_t i = 0; i < out.data.size(); ++i) out.data[i] = {re.data[i], im.data[i]};
    return out;
}

static RealImage2D real_part(const ComplexImage2D& c) {
    RealImage2D out = make_real(c.width, c.height);
    for (size_t i = 0; i < out.data.size(); ++i) out.data[i] = c.data[i].real();
    return out;
}
static RealImage2D imag_part(const ComplexImage2D& c) {
    RealImage2D out = make_real(c.width, c.height);
    for (size_t i = 0; i < out.data.size(); ++i) out.data[i] = c.data[i].imag();
    return out;
}

static RealImage2D add(const RealImage2D& a, const RealImage2D& b) {
    if (a.width != b.width || a.height != b.height) throw std::runtime_error("add size mismatch");
    RealImage2D out = make_real(a.width, a.height);
    for (size_t i = 0; i < out.data.size(); ++i) out.data[i] = a.data[i] + b.data[i];
    return out;
}
static RealImage2D sub(const RealImage2D& a, const RealImage2D& b) {
    if (a.width != b.width || a.height != b.height) throw std::runtime_error("sub size mismatch");
    RealImage2D out = make_real(a.width, a.height);
    for (size_t i = 0; i < out.data.size(); ++i) out.data[i] = a.data[i] - b.data[i];
    return out;
}
static RealImage2D scale(const RealImage2D& a, double s) {
    RealImage2D out = make_real(a.width, a.height);
    for (size_t i = 0; i < out.data.size(); ++i) out.data[i] = s * a.data[i];
    return out;
}

// Fuse two complex bands using local energy weighting (3x3 window)
// a: IR band, b: VIS band
static ComplexImage2D fuse_band_local_energy(const ComplexImage2D& a, const ComplexImage2D& b) {
    if (a.width != b.width || a.height != b.height) {
        throw std::runtime_error("fuse_band_local_energy: size mismatch");
    }
    const int w = a.width;
    const int h = a.height;
    const int R = 1;
    const double eps = 1e-8;

    ComplexImage2D out{w, h, {}};
    out.data.resize(a.data.size());

    auto mag = [](const std::complex<double>& c) { return std::abs(c); };

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double Ey = 0.0;
            double Ei = 0.0;
            for (int dy = -R; dy <= R; ++dy) {
                int yy = reflect_index(y + dy, h);
                for (int dx = -R; dx <= R; ++dx) {
                    int xx = reflect_index(x + dx, w);
                    const size_t idx_n = static_cast<size_t>(yy) * static_cast<size_t>(w) + static_cast<size_t>(xx);
                    Ey += std::fabs(mag(b.data[idx_n]));
                    Ei += std::fabs(mag(a.data[idx_n]));
                }
            }
            const double wi = Ei / (Ei + Ey + eps); // weight for IR (a)
            const double wy = 1.0 - wi;             // weight for VIS (b)
            const size_t idx = static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x);
            out.data[idx] = wy * b.data[idx] + wi * a.data[idx];
        }
    }
    return out;
}

// Fuse two real lowpass bands using local energy weighting (3x3 window)
// a: IR lowpass, b: VIS lowpass
static RealImage2D fuse_lowpass_local_energy(const RealImage2D& a, const RealImage2D& b) {
    if (a.width != b.width || a.height != b.height) {
        throw std::runtime_error("fuse_lowpass_local_energy: size mismatch");
    }
    const int w = a.width;
    const int h = a.height;
    const int R = 1;
    const double eps = 1e-8;

    RealImage2D out = make_real(w, h);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double Ey = 0.0;
            double Ei = 0.0;
            for (int dy = -R; dy <= R; ++dy) {
                int yy = reflect_index(y + dy, h);
                for (int dx = -R; dx <= R; ++dx) {
                    int xx = reflect_index(x + dx, w);
                    const size_t idx_n = static_cast<size_t>(yy) * static_cast<size_t>(w) + static_cast<size_t>(xx);
                    Ey += std::fabs(b.data[idx_n]);
                    Ei += std::fabs(a.data[idx_n]);
                }
            }
            const double wi = Ei / (Ei + Ey + eps); // weight for IR (a)
            const double wy = 1.0 - wi;             // weight for VIS (b)
            const size_t idx = static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x);
            out.data[idx] = wy * b.data[idx] + wi * a.data[idx];
        }
    }
    return out;
}

} // namespace

DtcwtPyramid dtcwt2d_forward(const RealImage2D& input, int levels)
{
    check_real_image(input);
    if (levels <= 0) throw std::runtime_error("dtcwt2d_forward: levels must be > 0");

    const FilterBank1D fb_near = near_sym_a_filterbank();
    const FilterBank1D fb_q = qshift10_filterbank();

    DtcwtPyramid pyr;
    pyr.levels.clear();
    pyr.levels.reserve(static_cast<size_t>(levels));

    // Dual-tree lowpass init
    RealImage2D lp_a = input;
    RealImage2D lp_b = input;

    const double s = 1.0 / std::sqrt(2.0);

    for (int j = 1; j <= levels; ++j) {
        // Select filters (level-dependent)
        // Per spec: level 1 uses near-symmetric, levels>=2 use Q-shift.
        const FilterBank1D& fb = (j == 1) ? fb_near : fb_q;

        // Build four trees AA, AB, BA, BB
        // Row: A shift=0, B shift=1 ; Col: A shift=0, B shift=1
        const Subbands2D aa = analysis_2d_tree(lp_a, fb.h0a, fb.h1a, 0, fb.h0a, fb.h1a, 0);
        const Subbands2D ab = analysis_2d_tree(lp_a, fb.h0a, fb.h1a, 0, fb.h0b, fb.h1b, 1);
        const Subbands2D ba = analysis_2d_tree(lp_b, fb.h0b, fb.h1b, 1, fb.h0a, fb.h1a, 0);
        const Subbands2D bb = analysis_2d_tree(lp_b, fb.h0b, fb.h1b, 1, fb.h0b, fb.h1b, 1);

        // Complex subband mixing (MANDATORY formulas; orientation order MUST match spec)
        // LH -> ±75°
        RealImage2D Re_p75 = scale(add(aa.LH, bb.LH), s);
        RealImage2D Im_p75 = scale(sub(ab.LH, ba.LH), s);
        RealImage2D Re_m75 = scale(sub(aa.LH, bb.LH), s);
        RealImage2D Im_m75 = scale(add(ab.LH, ba.LH), s);

        // HL -> ±15°
        RealImage2D Re_p15 = scale(add(aa.HL, bb.HL), s);
        RealImage2D Im_p15 = scale(sub(ba.HL, ab.HL), s);
        RealImage2D Re_m15 = scale(sub(aa.HL, bb.HL), s);
        RealImage2D Im_m15 = scale(add(ba.HL, ab.HL), s);

        // HH -> ±45°
        RealImage2D Re_p45 = scale(sub(aa.HH, bb.HH), s);
        RealImage2D Im_p45 = scale(add(ab.HH, ba.HH), s);
        RealImage2D Re_m45 = scale(add(aa.HH, bb.HH), s);
        RealImage2D Im_m45 = scale(sub(ba.HH, ab.HH), s);

        DtcwtLevel level;
        // Store per-level pre-decimation size (needed for exact inverse sizing)
        level.in_width = lp_a.width;
        level.in_height = lp_a.height;
        level.band[0] = pack_complex(Re_p15, Im_p15); // +15
        level.band[1] = pack_complex(Re_m15, Im_m15); // -15
        level.band[2] = pack_complex(Re_p45, Im_p45); // +45
        level.band[3] = pack_complex(Re_m45, Im_m45); // -45
        level.band[4] = pack_complex(Re_p75, Im_p75); // +75
        level.band[5] = pack_complex(Re_m75, Im_m75); // -75
        pyr.levels.push_back(std::move(level));

        // Lowpass propagation (CRITICAL): do NOT average
        lp_a = aa.LL; // Next lowpass_a = LL_AA
        lp_b = bb.LL; // Next lowpass_b = LL_BB
    }

    pyr.lowpass_a = std::move(lp_a);
    pyr.lowpass_b = std::move(lp_b);
    return pyr;
}

RealImage2D dtcwt2d_inverse(const DtcwtPyramid& pyramid)
{
    if (pyramid.levels.empty()) throw std::runtime_error("dtcwt2d_inverse: empty pyramid");
    check_real_image(pyramid.lowpass_a);
    check_real_image(pyramid.lowpass_b);

    // Start from coarsest lowpasses
    RealImage2D lp_a = pyramid.lowpass_a;
    RealImage2D lp_b = pyramid.lowpass_b;

    const FilterBank1D fb_near = near_sym_a_filterbank();
    const FilterBank1D fb_q = qshift10_filterbank();
    const double s = 1.0 / std::sqrt(2.0);
    const double inv = 1.0 / (2.0 * s); // per spec

    // Reconstruct from last level down to 1
    for (int j = static_cast<int>(pyramid.levels.size()); j >= 1; --j) {
        const DtcwtLevel& lvl = pyramid.levels[static_cast<size_t>(j - 1)];

        // Extract oriented bands
        const RealImage2D Re_p15 = real_part(lvl.band[0]);
        const RealImage2D Im_p15 = imag_part(lvl.band[0]);
        const RealImage2D Re_m15 = real_part(lvl.band[1]);
        const RealImage2D Im_m15 = imag_part(lvl.band[1]);

        const RealImage2D Re_p45 = real_part(lvl.band[2]);
        const RealImage2D Im_p45 = imag_part(lvl.band[2]);
        const RealImage2D Re_m45 = real_part(lvl.band[3]);
        const RealImage2D Im_m45 = imag_part(lvl.band[3]);

        const RealImage2D Re_p75 = real_part(lvl.band[4]);
        const RealImage2D Im_p75 = imag_part(lvl.band[4]);
        const RealImage2D Re_m75 = real_part(lvl.band[5]);
        const RealImage2D Im_m75 = imag_part(lvl.band[5]);

        // Inverse subband unmixing (EXACT inverse patterns, per spec)
        // From ±75° -> LH
        RealImage2D LH_AA = scale(add(Re_p75, Re_m75), inv);
        RealImage2D LH_BB = scale(sub(Re_p75, Re_m75), inv);
        RealImage2D LH_AB = scale(add(Im_p75, Im_m75), inv);
        RealImage2D LH_BA = scale(sub(Im_m75, Im_p75), inv);

        // From ±15° -> HL
        RealImage2D HL_AA = scale(add(Re_p15, Re_m15), inv);
        RealImage2D HL_BB = scale(sub(Re_p15, Re_m15), inv);
        RealImage2D HL_BA = scale(add(Im_p15, Im_m15), inv);
        RealImage2D HL_AB = scale(sub(Im_m15, Im_p15), inv);

        // From ±45° -> HH
        RealImage2D HH_AA = scale(add(Re_p45, Re_m45), inv);
        RealImage2D HH_BB = scale(sub(Re_m45, Re_p45), inv);
        RealImage2D HH_BA = scale(add(Im_p45, Im_m45), inv);
        RealImage2D HH_AB = scale(sub(Im_p45, Im_m45), inv);

        // Filters (level-dependent)
        const FilterBank1D& fb = (j == 1) ? fb_near : fb_q;

        // Need exact output size at this level (odd/even safe)
        if (lvl.in_width <= 0 || lvl.in_height <= 0) {
            throw std::runtime_error("dtcwt2d_inverse: missing per-level size metadata");
        }
        const int outW = lvl.in_width;
        const int outH = lvl.in_height;

        // Reconstruction
        RealImage2D recon_AA;
        RealImage2D recon_BB;
        if (j == 1) {
            // near-symmetric
            // Tree AA reconstruction (row A shift=0, col A shift=0)
            recon_AA = synth_2d_tree(lp_a, LH_AA, HL_AA, HH_AA,
                                                outW, outH,
                                                fb.g0a, fb.g1a, 0,
                                                fb.g0a, fb.g1a, 0);

            // Tree BB reconstruction (row B shift=1, col B shift=1)
            recon_BB = synth_2d_tree(lp_b, LH_BB, HL_BB, HH_BB,
                                                outW, outH,
                                                fb.g0b, fb.g1b, 1,
                                                fb.g0b, fb.g1b, 1);
        }else {
            // q-shift
            // Tree AA reconstruction (row A shift=0, col A shift=0)
            recon_AA = synth_2d_tree(lp_a, LH_AA, HL_AA, HH_AA,
                                                outW, outH,
                                                fb.g0a, fb.g1a, 0,   // row shift
                                                fb.g0a, fb.g1a, 0);  // col shift

            // Tree BB reconstruction (row B shift=1, col B shift=1)
            recon_BB = synth_2d_tree(lp_b, LH_BB, HL_BB, HH_BB,
                                                outW, outH,
                                                fb.g0b, fb.g1b, 1,   // row shift
                                                fb.g0b, fb.g1b, 1);  // col shift
        }

        // CRITICAL (per spec): AA/BB recon are intermediate.
        // True lowpass to propagate to next (finer) stage is:
        // lp = 0.5 * (recon_AA + recon_BB)
        RealImage2D lp = make_real(outW, outH);
        for (size_t i = 0; i < lp.data.size(); ++i) {
            lp.data[i] = 0.5 * (recon_AA.data[i] + recon_BB.data[i]);
        }

        // Prepare next lowpass for the next (finer) stage:
        // both trees use the same propagated lowpass (do NOT keep them split).
        lp_a = lp;
        lp_b = lp;

        if (j == 1) {
            return lp;
        }
    }

    throw std::runtime_error("dtcwt2d_inverse: unreachable");
}

bool dtcwt_fused_image(const std::string& vis_path, 
    const std::string& ir_path,
    const std::string& output_path)
{
    const int levels = 1;

    // 讀圖：VIS RGB，IR 灰階
    Image vis_img;
    Image ir_img;
    if (!load_image(vis_path, 3, vis_img)) {
        std::cerr << "[Error] Cannot load VIS image: " << vis_path << "\n";
        return false;
    }
    if (!load_image(ir_path, 1, ir_img)) {
        std::cerr << "[Error] Cannot load IR image: " << ir_path << "\n";
        return false;
    }

    if (vis_img.width != ir_img.width || vis_img.height != ir_img.height) {
        std::cerr << "[Error] VIS/IR size mismatch! VIS="
                  << vis_img.width << "x" << vis_img.height
                  << ", IR=" << ir_img.width << "x" << ir_img.height << "\n";
        return false;
    }

    // RGB -> YCbCr；IR -> double
    std::vector<double> Y, Cb, Cr, IR;
    rgb_to_ycbcr(vis_img, Y, Cb, Cr);
    ir_to_double(ir_img, IR);

    // 準備 RealImage2D（Y、IR 已是 double）
    RealImage2D vis{vis_img.width, vis_img.height, std::move(Y)};
    RealImage2D ir{ir_img.width, ir_img.height, std::move(IR)};

    // forward
    DtcwtPyramid pyr_vis = dtcwt2d_forward(vis, levels);
    DtcwtPyramid pyr_ir  = dtcwt2d_forward(ir,  levels);

    // fuse
    DtcwtPyramid pyr_fused;
    pyr_fused.levels.resize(pyr_ir.levels.size());
    pyr_fused.lowpass_a = fuse_lowpass_local_energy(pyr_ir.lowpass_a, pyr_vis.lowpass_a);
    pyr_fused.lowpass_b = fuse_lowpass_local_energy(pyr_ir.lowpass_b, pyr_vis.lowpass_b);
    for (size_t lev = 0; lev < pyr_ir.levels.size(); ++lev) {
        // 保留尺寸 metadata，供 inverse 使用
        pyr_fused.levels[lev].in_width = pyr_ir.levels[lev].in_width;
        pyr_fused.levels[lev].in_height = pyr_ir.levels[lev].in_height;
        for (int k = 0; k < 6; ++k) {
            const ComplexImage2D& a = pyr_ir.levels[lev].band[k];
            const ComplexImage2D& b = pyr_vis.levels[lev].band[k];
            if (a.width != b.width || a.height != b.height) {
                std::cerr << "[Error] Subband size mismatch at level=" << lev << " band=" << k
                            << " IR=" << a.width << "x" << a.height
                            << " VIS=" << b.width << "x" << b.height << "\n";
                return 1;
            }
            pyr_fused.levels[lev].band[k] = fuse_band_local_energy(a, b);
        }
    }

    // inverse
    RealImage2D fused = dtcwt2d_inverse(pyr_fused);

    // Y_fused + 原本 Cb/Cr -> fused RGB
    Image fused_img;
    fused_img.width = vis_img.width;
    fused_img.height = vis_img.height;
    ycbcr_to_rgb(fused.data, Cb, Cr, fused_img);

    if (!write_image_png(output_path, fused_img)) {
        std::cerr << "[Error] Failed to write: " << output_path << "\n";
        return false;
    }

    return true;
}


int run_dtcwt_single(const Options& opt) {
    std::cout << "[Info] Mode: dtcwt (single)\n";
    std::cout << "[Info] VIS: " << opt.vis_path << "\n";
    std::cout << "[Info] IR : " << opt.ir_path << "\n";
    std::cout << "[Info] Output directory: " << opt.output_dir << "\n";

    if (!fs::exists(opt.output_dir)) {
        fs::create_directories(opt.output_dir);
    }
    const std::string filename = fs::path(opt.vis_path).stem().string();
    const std::string output_path = opt.output_dir + "/fused_" + filename + "_" + opt.mode + ".png";

    bool ok = dtcwt_fused_image(opt.vis_path, opt.ir_path, output_path);

    if (ok) std::cout << "[Info] Success\n";
    return ok ? 0 : 1; 
}

int run_dtcwt_batch(const Options& opt) {
    std::cout << "[Info] Mode: dtcwt (batch)\n";
    std::cout << "[Info] VIS directory: " << opt.vis_dir << "\n";
    std::cout << "[Info] IR directory: " << opt.ir_dir << "\n";
    std::cout << "[Info] Output directory: " << opt.output_dir << "\n";

    if (!fs::exists(opt.output_dir)) {
        fs::create_directories(opt.output_dir);
    }

    std::vector<std::string> vis_files = get_png_files(opt.vis_dir);
    std::vector<std::string> ir_files = get_png_files(opt.ir_dir);
    if (vis_files.empty()) {
        std::cerr << "[Error] No PNG images found in VIS directory\n";
        return 1;
    }
    if (ir_files.empty()) {
        std::cerr << "[Error] No PNG images found in IR directory\n";
        return 1;
    }

    std::unordered_set<std::string> ir_set;
    ir_set.reserve(ir_files.size() * 2);
    for (const auto& s : ir_files) ir_set.insert(s);

    int success_count = 0;
    int fail_count = 0;

    for (const auto& filename : vis_files) {
        if (!ir_set.count(filename)) {
            std::cout << "No corresponding IR image: " << filename << ".png, skip\n";
            fail_count++;
            continue;
        }

        std::string vis_path = opt.vis_dir + "/" + filename + ".png";
        std::string ir_path = opt.ir_dir + "/" + filename + ".png";
        std::string output_path = opt.output_dir + "/fused_" + filename + "_" + opt.mode + ".png";

        std::cout << "[" << (success_count + fail_count + 1) << "/" << vis_files.size()
            << "] Processing: " << filename << ".png ... ";

        if (dtcwt_fused_image(vis_path, ir_path, output_path)) {
            std::cout << "Success\n";
            success_count++;
        } else {
            std::cout << "Failed\n";
            fail_count++;
        }
    }

    std::cout << "========================================\n";
    std::cout << "[Info] Completed! Success: " << success_count
        << ", Failed: " << fail_count << "\n";

    return (fail_count > 0) ? 1 : 0;
}


int test_dtcwt_single(const Options& opt) {
    try {
        const int levels = 1; // 可依需求調整
        std::cout << "[Info] Mode: TEST dtcwt (single)\n";
        std::cout << "[Info] VIS: " << opt.vis_path << "\n";
        std::cout << "[Info] IR : " << opt.ir_path << "\n";
        std::cout << "[Info] Output directory: " << opt.output_dir << "\n";
        std::cout << "[Info] Levels: " << levels << "\n";

        if (!fs::exists(opt.output_dir)) {
            fs::create_directories(opt.output_dir);
        }
        const std::string filename = fs::path(opt.vis_path).stem().string();
        const std::string output_path = opt.output_dir + "/fused_" + filename + "_" + opt.mode + ".png";

        // 讀圖：VIS RGB，IR 灰階
        Image vis_img;
        Image ir_img;
        if (!load_image(opt.vis_path, 3, vis_img)) {
            std::cerr << "[Error] Cannot load VIS image: " << opt.vis_path << "\n";
            return 1;
        }
        if (!load_image(opt.ir_path, 1, ir_img)) {
            std::cerr << "[Error] Cannot load IR image: " << opt.ir_path << "\n";
            return 1;
        }

        if (vis_img.width != ir_img.width || vis_img.height != ir_img.height) {
            std::cerr << "[Error] VIS/IR size mismatch! VIS="
                      << vis_img.width << "x" << vis_img.height
                      << ", IR=" << ir_img.width << "x" << ir_img.height << "\n";
            return 1;
        }

        // RGB -> YCbCr；IR -> double
        std::vector<double> Y, Cb, Cr, IR;
        rgb_to_ycbcr(vis_img, Y, Cb, Cr);
        ir_to_double(ir_img, IR);

        // 準備 RealImage2D（Y、IR 已是 double）
        RealImage2D vis{vis_img.width, vis_img.height, std::move(Y)};
        RealImage2D ir{ir_img.width, ir_img.height, std::move(IR)};

        // quick stats helper
        auto print_stats = [](const char* tag, const RealImage2D& im) {
            double mn = im.data.empty() ? 0.0 : im.data[0];
            double mx = im.data.empty() ? 0.0 : im.data[0];
            long double sum = 0.0L;
            for (double v : im.data) {
                if (v < mn) mn = v;
                if (v > mx) mx = v;
                sum += static_cast<long double>(v);
            }
            const double mean = im.data.empty() ? 0.0 : static_cast<double>(sum / static_cast<long double>(im.data.size()));
            std::cerr << "[DTCWT][STATS] " << tag << " min=" << mn << " max=" << mx << " mean=" << mean << "\n";
        };
        print_stats("IR input", ir);
        print_stats("VIS input (Y)", vis);

        // forward
        std::cerr << "[DTCWT] forward VIS levels=" << levels << "\n";
        DtcwtPyramid pyr_vis = dtcwt2d_forward(vis, levels);
        std::cerr << "[DTCWT] forward IR levels=" << levels << "\n";
        DtcwtPyramid pyr_ir  = dtcwt2d_forward(ir,  levels);  fs::create_directories(opt.output_dir);

         // inverse
        std::cerr << "[DTCWT] inverse\n";
        RealImage2D ir_recon = dtcwt2d_inverse(pyr_ir);
        print_stats("IR inverse", ir_recon);
        RealImage2D vis_recon = dtcwt2d_inverse(pyr_vis);
        print_stats("VIS inverse (Y)", vis_recon);
        
        // 將 IR 和 VIS (Y) 寫出為圖片
        // IR: 單通道灰階
        Image ir_out;
        ir_out.width = ir_recon.width;
        ir_out.height = ir_recon.height;
        ir_out.channels = 1;
        ir_out.data.resize(ir_recon.data.size());

        for (size_t i = 0; i < ir_recon.data.size(); ++i) {
            double val = std::clamp(ir_recon.data[i], 0.0, 255.0);
            ir_out.data[i] = static_cast<unsigned char>(std::round(val));
        }
        std::string ir_output_path = opt.output_dir + "/ir_inverse_" + filename + ".png";
        if (!write_image_png(ir_output_path, ir_out)) {
            std::cerr << "[Error] Failed to write IR output: " << ir_output_path << "\n";
            return 1;
        }
        std::cout << "[Info] IR output written to: " << ir_output_path << "\n";

        // VIS: Y 通道轉回 RGB (保持 Cb, Cr 不變)
        Image vis_out;
        vis_out.width = vis_recon.width;
        vis_out.height = vis_recon.height;
        vis_out.channels = 3;
        vis_out.data.resize(vis_recon.width * vis_recon.height * 3);

        ycbcr_to_rgb(vis_recon.data, Cb, Cr, vis_out);
        std::string vis_output_path = opt.output_dir + "/vis_inverse_" + filename + ".png";
        if (!write_image_png(vis_output_path, vis_out)) {
            std::cerr << "[Error] Failed to write VIS output: " << vis_output_path << "\n";
            return 1;
        }
        std::cout << "[Info] VIS output written to: " << vis_output_path << "\n";
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[DTCWT][EXCEPTION] " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "[DTCWT][EXCEPTION] unknown exception\n";
        return 1;
    }
}
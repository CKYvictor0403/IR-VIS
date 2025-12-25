#pragma once

#include <array>
#include <complex>
#include <vector>
#include "utils.h"

struct RealImage2D {
    int width;
    int height;
    std::vector<double> data; // row-major
};

struct ComplexImage2D {
    int width;
    int height;
    std::vector<std::complex<double>> data;
};

struct DtcwtLevel {
    // Size of the (per-level) input lowpass BEFORE decimation at this level.
    // This is REQUIRED so inverse can reconstruct exact odd/even dimensions.
    int in_width = 0;
    int in_height = 0;

    // Orientation order MUST be preserved
    // [0] +15°, [1] -15°, [2] +45°, [3] -45°, [4] +75°, [5] -75°
    std::array<ComplexImage2D, 6> band;
};

struct DtcwtPyramid {
    std::vector<DtcwtLevel> levels;
    // Dual-tree lowpass is REQUIRED for invertibility
    RealImage2D lowpass_a;
    RealImage2D lowpass_b;
};

DtcwtPyramid dtcwt2d_forward(const RealImage2D& input, int levels);
RealImage2D dtcwt2d_inverse(const DtcwtPyramid& pyramid);

// PR self-test (single image, no file output): forward -> inverse -> print stats
int test_dtcwt_single(const Options& opt);

int run_dtcwt_single(const Options& opt);
int run_dtcwt_batch(const Options& opt);
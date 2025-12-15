#include "dwt.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include "wavelib.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace fs = std::filesystem;

struct Image {
    int width = 0;
    int height = 0;
    int channels = 0; // e.g. 3 = RGB, 1 = Gray
    std::vector<unsigned char> data;
};

static bool load_image(const std::string& path,
    int desired_channels,
    Image& out_img)
{
    int w, h, c;
    stbi_uc* pixels = stbi_load(path.c_str(), &w, &h, &c, desired_channels);
    if (!pixels) {
        std::cerr << "[Error] Failed to load image: " << path << "\n"
            << "        stbi failure reason: " << stbi_failure_reason() << std::endl;
        return false;
    }

    out_img.width = w;
    out_img.height = h;
    out_img.channels = (desired_channels > 0) ? desired_channels : c;
    out_img.data.assign(pixels, pixels + w * h * out_img.channels);

    stbi_image_free(pixels);
    return true;
}

std::vector<std::string> get_png_files(const std::string& dir_path)
{
    std::vector<std::string> files;
    if (!fs::exists(dir_path) || !fs::is_directory(dir_path)) {
        std::cerr << "[Error] Directory does not exist: " << dir_path << "\n";
        return files;
    }

    for (const auto& entry : fs::directory_iterator(dir_path)) {
        if (entry.is_regular_file()) {
            std::string ext = entry.path().extension().string();
            std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
            if (ext == ".png") {
                files.push_back(entry.path().filename().stem().string());
            }
        }
    }

    std::sort(files.begin(), files.end());
    return files;
}

bool get_image_dimensions(const std::string& path, int& out_width, int& out_height)
{
    int channels = 0;
    if (!stbi_info(path.c_str(), &out_width, &out_height, &channels)) {
        return false;
    }
    return true;
}

static bool write_image_png(const std::string& path, const Image& img)
{
    if (img.data.empty()) {
        std::cerr << "[Error] write_image_png: empty image data\n";
        return false;
    }

    int stride_in_bytes = img.width * img.channels;
    int ok = stbi_write_png(path.c_str(),
        img.width,
        img.height,
        img.channels,
        img.data.data(),
        stride_in_bytes);

    if (!ok) {
        std::cerr << "[Error] Failed to write PNG: " << path << std::endl;
        return false;
    }
    return true;
}

static void rgb_to_ycbcr(const Image& rgb,
    std::vector<double>& Y,
    std::vector<double>& Cb,
    std::vector<double>& Cr)
{
    const int w = rgb.width;
    const int h = rgb.height;
    Y.resize(w * h);
    Cb.resize(w * h);
    Cr.resize(w * h);

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int idx = (y * w + x) * 3;
            double R = rgb.data[idx + 0];
            double G = rgb.data[idx + 1];
            double B = rgb.data[idx + 2];

            double Yv = 0.299 * R + 0.587 * G + 0.114 * B;
            double Cbv = (B - Yv) / 1.772 + 128.0;
            double Crv = (R - Yv) / 1.402 + 128.0;

            Y[y * w + x] = Yv;
            Cb[y * w + x] = Cbv;
            Cr[y * w + x] = Crv;
        }
    }
}

static void ycbcr_to_rgb(const std::vector<double>& Y,
    const std::vector<double>& Cb,
    const std::vector<double>& Cr,
    Image& rgb)
{
    const int w = rgb.width;
    const int h = rgb.height;
    rgb.channels = 3;
    rgb.data.resize(w * h * 3);

    auto clamp255 = [](double v) -> unsigned char {
        if (v < 0.0) v = 0.0;
        if (v > 255.0) v = 255.0;
        return static_cast<unsigned char>(v + 0.5);
        };

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int idx = y * w + x;
            double Yv = Y[idx];
            double Cbv = Cb[idx];
            double Crv = Cr[idx];

            double R = Yv + 1.402 * (Crv - 128.0);
            double G = Yv - 0.344136 * (Cbv - 128.0)
                - 0.714136 * (Crv - 128.0);
            double B = Yv + 1.772 * (Cbv - 128.0);

            int idx3 = idx * 3;
            rgb.data[idx3 + 0] = clamp255(R);
            rgb.data[idx3 + 1] = clamp255(G);
            rgb.data[idx3 + 2] = clamp255(B);
        }
    }
}

static void ir_to_double(const Image& ir, std::vector<double>& out)
{
    const int w = ir.width;
    const int h = ir.height;
    out.resize(w * h);
    for (int i = 0; i < w * h; ++i) {
        out[i] = static_cast<double>(ir.data[i]);
    }
}

static void fuse_subband_local_energy(const double* Yc,
    const double* IRc,
    double* F,
    int rows,
    int cols,
    bool is_LL,
    double enhancement_factor = 1.0)
{
    const int R = 1;
    const double eps = 1e-8;

    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            double Ey = 0.0;
            double Ei = 0.0;

            for (int dr = -R; dr <= R; ++dr) {
                int rr = r + dr;
                if (rr < 0 || rr >= rows) continue;
                for (int dc = -R; dc <= R; ++dc) {
                    int cc = c + dc;
                    if (cc < 0 || cc >= cols) continue;

                    double vy = Yc[rr * cols + cc];
                    double vi = IRc[rr * cols + cc];
                    Ey += fabs(vy);
                    Ei += fabs(vi);
                }
            }

            double wy, wi;
            if (is_LL) {
                wi = Ei / (Ei + Ey + eps);
                wy = 1.0 - wi;
            } else {
                wi = Ei / (Ei + Ey + eps);
                wy = 1.0 - wi;
            }

            int idx = r * cols + c;
            F[idx] = wy * Yc[idx] + wi * IRc[idx];

            if (!is_LL) {
                F[idx] *= enhancement_factor;
            }
        }
    }
}

// 融合單張圖片的函式
bool dwt_fused_image(const std::string& rgb_path, 
    const std::string& ir_path,
    const std::string& output_path,
    wt2_object wt)
{
    Image rgb_img;
    Image ir_img;

    // 讀取 RGB：期望 3 channel（R,G,B）
    if (!load_image(rgb_path, 3, rgb_img)) {
        std::cerr << "[Error] Cannot load RGB image: " << rgb_path << "\n";
        return false;
    }

    // 讀取 IR：通常是 1 channel 灰階
    if (!load_image(ir_path, 1, ir_img)) {
        std::cerr << "[Error] Cannot load IR image: " << ir_path << "\n";
        return false;
    }

    // 檢查兩張圖尺寸是否一致
    if (rgb_img.width != ir_img.width || rgb_img.height != ir_img.height) {
        std::cerr << "[Error] RGB and IR size mismatch! RGB: " 
        << rgb_img.width << "x" << rgb_img.height 
        << ", IR: " << ir_img.width << "x" << ir_img.height << "\n";
        return false;
    }

    const int w = rgb_img.width;
    const int h = rgb_img.height;

    // 檢查 DWT 物件尺寸是否與輸入一致（wt2_object 是固定 rows/cols 的）
    if (wt == nullptr) {
        std::cerr << "[Error] wt2_object is null\n";
        return false;
    }
    if (wt->rows != h || wt->cols != w) {
        std::cerr << "[Error] wt2_object size mismatch! wt=(" << wt->cols << "x" << wt->rows
        << "), img=(" << w << "x" << h << ")\n";
        return false;
    }

    // 1) RGB -> YCbCr、IR -> double
    std::vector<double> Y, Cb, Cr, IR;
    rgb_to_ycbcr(rgb_img, Y, Cb, Cr);
    ir_to_double(ir_img, IR);

    // 3) 對 Y 與 IR 各做一次 2D DWT
    double* coeff_Y = dwt2(wt, Y.data());
    double* coeff_IR = dwt2(wt, IR.data());

    // 4) 準備一個 fused 的係數向量
    std::vector<double> coeff_fused(wt->outlength, 0.0);

    int sub_rows = 0, sub_cols = 0;

    // --- LL (Approximation, "A")：以 IR 為主，局部能量加權 ---
    {
        double* LL_Y = getWT2Coeffs(wt, coeff_Y, 1, (char*)"A", &sub_rows, &sub_cols);
        double* LL_IR = getWT2Coeffs(wt, coeff_IR, 1, (char*)"A", &sub_rows, &sub_cols);
        double* LL_F = getWT2Coeffs(wt, coeff_fused.data(), 1, (char*)"A", &sub_rows, &sub_cols);

        fuse_subband_local_energy(LL_Y, LL_IR, LL_F, sub_rows, sub_cols, /*is_LL=*/true);
    }

    // --- LH (Horizontal detail, "H") ---
    {
        double* LH_Y = getWT2Coeffs(wt, coeff_Y, 1, (char*)"H", &sub_rows, &sub_cols);
        double* LH_IR = getWT2Coeffs(wt, coeff_IR, 1, (char*)"H", &sub_rows, &sub_cols);
        double* LH_F = getWT2Coeffs(wt, coeff_fused.data(), 1, (char*)"H", &sub_rows, &sub_cols);

        fuse_subband_local_energy(LH_Y, LH_IR, LH_F, sub_rows, sub_cols, /*is_LL=*/false, 1.25);
    }

    // --- HL (Vertical detail, "V") ---
    {
        double* HL_Y = getWT2Coeffs(wt, coeff_Y, 1, (char*)"V", &sub_rows, &sub_cols);
        double* HL_IR = getWT2Coeffs(wt, coeff_IR, 1, (char*)"V", &sub_rows, &sub_cols);
        double* HL_F = getWT2Coeffs(wt, coeff_fused.data(), 1, (char*)"V", &sub_rows, &sub_cols);

        fuse_subband_local_energy(HL_Y, HL_IR, HL_F, sub_rows, sub_cols, /*is_LL=*/false, 1.25);
    }

    // --- HH (Diagonal detail, "D") ---
    {
        double* HH_Y = getWT2Coeffs(wt, coeff_Y, 1, (char*)"D", &sub_rows, &sub_cols);
        double* HH_IR = getWT2Coeffs(wt, coeff_IR, 1, (char*)"D", &sub_rows, &sub_cols);
        double* HH_F = getWT2Coeffs(wt, coeff_fused.data(), 1, (char*)"D", &sub_rows, &sub_cols);

        fuse_subband_local_energy(HH_Y, HH_IR, HH_F, sub_rows, sub_cols, /*is_LL=*/false, 1.5);
    }

    // 5) 反變換，重建 Y_fused
    std::vector<double> Y_fused(w * h);
    idwt2(wt, coeff_fused.data(), Y_fused.data());

    // 清理 DWT 資源
    free(coeff_Y);
    free(coeff_IR);

    // 6) Y_fused + 原本 Cb/Cr -> fused RGB
    Image fused_img;
    fused_img.width = w;
    fused_img.height = h;

    ycbcr_to_rgb(Y_fused, Cb, Cr, fused_img);

    // 7) 寫出融合結果
    if (!write_image_png(output_path, fused_img)) {
        std::cerr << "[Error] Failed to write: " << output_path << "\n";
        return false;
    }

    return true;
}

int run_dwt_single(const Options& opt) {
    const std::string wave_name = "sym4";
    std::cout << "[Info] Mode: dwt (single)\n";
    std::cout << "[Info] VIS: " << opt.vis_path << "\n";
    std::cout << "[Info] IR : " << opt.ir_path << "\n";
    std::cout << "[Info] Output directory: " << opt.output_dir << "\n";
    std::cout << "[Info] Wave: " << wave_name << "\n";

    int w = 0, h = 0;
    if (!get_image_dimensions(opt.vis_path, w, h)) {
        std::cerr << "[Error] Cannot read VIS image size: " << opt.vis_path << "\n";
        return 1;
    }

    if (!fs::exists(opt.output_dir)) {
        fs::create_directories(opt.output_dir);
    }
    std::string filename = fs::path(opt.vis_path).stem().string();
    std::string output_path = opt.output_dir + "/fused_" + filename + "_" + opt.mode + ".png";

    wave_object wave = wave_init(const_cast<char*>(wave_name.c_str()));
    if (!wave) {
        std::cerr << "[Error] wave_init failed\n";
        return 1;
    }
    wt2_object wt = wt2_init(wave, const_cast<char*>("dwt"), h, w, 1);
    if (!wt) {
        std::cerr << "[Error] wt2_init failed\n";
        wave_free(wave);
        return 1;
    }

    bool ok = dwt_fused_image(opt.vis_path, opt.ir_path, output_path, wt);
    if (ok) std::cout << "[Info] Success\n";

    wt2_free(wt);
    wave_free(wave);
    return ok ? 0 : 1;
}

int run_dwt_batch(const Options& opt) {
    const std::string wave_name = "sym4";
    std::cout << "[Info] Mode: dwt (batch)\n";
    std::cout << "[Info] VIS directory: " << opt.vis_dir << "\n";
    std::cout << "[Info] IR directory: " << opt.ir_dir << "\n";
    std::cout << "[Info] Output directory: " << opt.output_dir << "\n";
    std::cout << "[Info] Wave: " << wave_name << "\n";

    if (!fs::exists(opt.output_dir)) {
        fs::create_directories(opt.output_dir);
        std::cout << "[Info] 建立輸出資料夾: " << opt.output_dir << "\n";
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

    wave_object wave = nullptr;
    wt2_object wt = nullptr;
    int current_w = 0;
    int current_h = 0;

    auto reset_dwt_context = [&]() {
        if (wt) {
            wt2_free(wt);
            wt = nullptr;
        }
        if (wave) {
            wave_free(wave);
            wave = nullptr;
        }
    };

    auto ensure_dwt_context = [&](int width, int height) -> bool {
        if (wt && wave && width == current_w && height == current_h) {
            return true;
        }
        reset_dwt_context();
        wave = wave_init(const_cast<char*>(wave_name.c_str()));
        if (!wave) {
            std::cerr << "[Error] Cannot initialize wave (" << wave_name << ")\n";
            return false;
        }
        wt = wt2_init(wave, const_cast<char*>("dwt"), height, width, 1);
        if (!wt) {
            std::cerr << "[Error] Cannot initialize wt2_object (" << width << "x" << height << ")\n";
            wave_free(wave);
            wave = nullptr;
            return false;
        }
        current_w = width;
        current_h = height;
        return true;
    };

    for (const auto& filename : vis_files) {
        if (!ir_set.count(filename)) {
            std::cout << "No corresponding IR image: " << filename << ".png, skip\n";
            fail_count++;
            continue;
        }

        std::string vis_path = opt.vis_dir + "/" + filename + ".png";
        std::string ir_path = opt.ir_dir + "/" + filename + ".png";
        std::string output_path = opt.output_dir + "/fused_" + filename + "_" + opt.mode + ".png";

        int vis_w = 0, vis_h = 0;
        if (!get_image_dimensions(vis_path, vis_w, vis_h)) {
            std::cerr << "[Warning] Cannot read VIS image size: " << vis_path << ", skip\n";
            fail_count++;
            continue;
        }
        if (!ensure_dwt_context(vis_w, vis_h)) {
            reset_dwt_context();
            return 1;
        }

        std::cout << "[" << (success_count + fail_count + 1) << "/" << vis_files.size()
            << "] Processing: " << filename << ".png ... ";

        if (dwt_fused_image(vis_path, ir_path, output_path, wt)) {
            std::cout << "Success\n";
            success_count++;
        } else {
            std::cout << "Failed\n";
            fail_count++;
        }
    }

    reset_dwt_context();

    std::cout << "========================================\n";
    std::cout << "[Info] Completed! Success: " << success_count
        << ", Failed: " << fail_count << "\n";

    return (fail_count > 0) ? 1 : 0;
}

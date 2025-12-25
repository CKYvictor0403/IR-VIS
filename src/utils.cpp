#include "utils.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>


#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

namespace fs = std::filesystem;

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

bool load_image(const std::string& path,
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

bool write_image_png(const std::string& path, const Image& img)
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

void rgb_to_ycbcr(const Image& rgb,
                  std::vector<double>& Y,
                  std::vector<double>& Cb,
                  std::vector<double>& Cr)
{
    const int w = rgb.width;
    const int h = rgb.height;
    Y.resize(static_cast<size_t>(w) * static_cast<size_t>(h));
    Cb.resize(static_cast<size_t>(w) * static_cast<size_t>(h));
    Cr.resize(static_cast<size_t>(w) * static_cast<size_t>(h));

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const int idx = (y * w + x) * 3;
            const double R = rgb.data[idx + 0];
            const double G = rgb.data[idx + 1];
            const double B = rgb.data[idx + 2];

            const double Yv = 0.299 * R + 0.587 * G + 0.114 * B;
            const double Cbv = (B - Yv) / 1.772 + 128.0;
            const double Crv = (R - Yv) / 1.402 + 128.0;

            Y[static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)] = Yv;
            Cb[static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)] = Cbv;
            Cr[static_cast<size_t>(y) * static_cast<size_t>(w) + static_cast<size_t>(x)] = Crv;
        }
    }
}

void ycbcr_to_rgb(const std::vector<double>& Y,
                  const std::vector<double>& Cb,
                  const std::vector<double>& Cr,
                  Image& rgb)
{
    const int w = rgb.width;
    const int h = rgb.height;
    rgb.channels = 3;
    rgb.data.resize(static_cast<size_t>(w) * static_cast<size_t>(h) * 3);

    auto clamp255 = [](double v) -> unsigned char {
        if (v < 0.0) v = 0.0;
        if (v > 255.0) v = 255.0;
        return static_cast<unsigned char>(std::round(v));
    };

    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            const int idx = y * w + x;
            const double Yv = Y[static_cast<size_t>(idx)];
            const double Cbv = Cb[static_cast<size_t>(idx)];
            const double Crv = Cr[static_cast<size_t>(idx)];

            const double R = Yv + 1.402 * (Crv - 128.0);
            const double G = Yv - 0.344136 * (Cbv - 128.0) - 0.714136 * (Crv - 128.0);
            const double B = Yv + 1.772 * (Cbv - 128.0);

            const int idx3 = idx * 3;
            rgb.data[static_cast<size_t>(idx3) + 0] = clamp255(R);
            rgb.data[static_cast<size_t>(idx3) + 1] = clamp255(G);
            rgb.data[static_cast<size_t>(idx3) + 2] = clamp255(B);
        }
    }
}

void ir_to_double(const Image& ir, std::vector<double>& out)
{
    const int w = ir.width;
    const int h = ir.height;
    out.resize(static_cast<size_t>(w) * static_cast<size_t>(h));
    for (int i = 0; i < w * h; ++i) {
        out[static_cast<size_t>(i)] = static_cast<double>(ir.data[static_cast<size_t>(i)]);
    }
}

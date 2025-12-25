#pragma once

#include <string>
#include <vector>

// 共用的命令列參數結構
struct Options {
    std::string mode;
    bool batch = false;
    std::string vis_dir;
    std::string ir_dir;
    std::string output_dir; // MUST be provided via --output <dir>
    bool output_provided = false;
    std::string vis_path;
    std::string ir_path;
};

// 影像基本結構
struct Image {
    int width = 0;
    int height = 0;
    int channels = 0; // e.g. 3 = RGB, 1 = Gray
    std::vector<unsigned char> data;
};

// Shared helpers
std::vector<std::string> get_png_files(const std::string& dir_path);
bool get_image_dimensions(const std::string& path, int& out_width, int& out_height);
bool load_image(const std::string& path, int desired_channels, Image& out_img);
bool write_image_png(const std::string& path, const Image& img);
void rgb_to_ycbcr(const Image& rgb, std::vector<double>& Y, std::vector<double>& Cb, std::vector<double>& Cr);
void ycbcr_to_rgb(const std::vector<double>& Y, const std::vector<double>& Cb, const std::vector<double>& Cr, Image& rgb);
void ir_to_double(const Image& ir, std::vector<double>& out);


#include <iostream>
#include <string>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace std;

// 簡單包一個 Image 結構，之後做 DWT / Fusion 比較好用
struct Image {
    int width = 0;
    int height = 0;
    int channels = 0; // e.g. 3 = RGB, 1 = Gray
    std::vector<unsigned char> data; // size = width * height * channels
};

// 讀圖工具函式
// desired_channels = 0 代表維持原通道數
// 例如 RGB 圖用 3、IR 灰階圖用 1
bool load_image(const std::string& path,
    int desired_channels,
    Image& out_img)
{
    int w, h, c;
    // stbi_uc 就是 unsigned char
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

    // 原始記憶體釋放
    stbi_image_free(pixels);
    return true;
}

// 簡單寫圖工具函式（方便 debug 用）
// 將 Image 寫成 PNG
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

int main(int argc, char** argv)
{
    // 你可以用命令列參數指定路徑：
    //   FinalProject.exe rgb.png ir.png
    // 若沒給，就用預設檔名
    std::string rgb_path = "C:/git_repos/MMIP/FinalProject/assets/vis/00000.png";
    std::string ir_path = "C:/git_repos/MMIP/FinalProject/assets/ir/00000.png";

    if (argc >= 3) {
        rgb_path = argv[1];
        ir_path = argv[2];
    }

    std::cout << "[Info] RGB image path: " << rgb_path << "\n";
    std::cout << "[Info] IR  image path: " << ir_path << "\n";

    Image rgb_img;
    Image ir_img;

    // 讀取 RGB：期望 3 channel（R,G,B）
    if (!load_image(rgb_path, 3, rgb_img)) {
        std::cerr << "[Fatal] Cannot load RGB image. Abort.\n";
        return 1;
    }

    // 讀取 IR：通常是 1 channel 灰階
    if (!load_image(ir_path, 1, ir_img)) {
        std::cerr << "[Fatal] Cannot load IR image. Abort.\n";
        return 1;
    }

    std::cout << "[Info] RGB  size: " << rgb_img.width << " x "
        << rgb_img.height << " x " << rgb_img.channels << "\n";
    std::cout << "[Info] IR   size: " << ir_img.width << " x "
        << ir_img.height << " x " << ir_img.channels << "\n";

    // ===== 尺寸檢查：是否都是 1024 x 768 =====
    const int EXPECT_W = 1024;
    const int EXPECT_H = 768;

    if (rgb_img.width != EXPECT_W || rgb_img.height != EXPECT_H) {
        std::cerr << "[Warning] RGB image is not 1024x768 ("
            << rgb_img.width << "x" << rgb_img.height << ")\n";
    }
    if (ir_img.width != EXPECT_W || ir_img.height != EXPECT_H) {
        std::cerr << "[Warning] IR image is not 1024x768 ("
            << ir_img.width << "x" << ir_img.height << ")\n";
    }

    // 再檢查兩張圖尺寸是否一致
    if (rgb_img.width != ir_img.width || rgb_img.height != ir_img.height) {
        std::cerr << "[Fatal] RGB and IR image sizes do not match.\n";
        return 1;
    }
    std::cout << "[Info] Images loaded successfully. Ready for DWT / fusion.\n";
    
    // =====  DWT fusion =====







    // 寫一份輸出確認讀取正常
    if (!write_image_png("build/output/rgb_check.png", rgb_img)) {
        std::cerr << "[Warning] Failed to write rgb_check.png\n";
    }
    if (!write_image_png("build/output/ir_check.png", ir_img)) {
        std::cerr << "[Warning] Failed to write ir_check.png\n";
    }

    return 0;
}
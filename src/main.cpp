#include <iostream>
#include <string>
#include <vector>

#include "dwt.h"

static bool parse_args(int argc, char** argv, Options& opt) {
    if (argc < 2) return false;

    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--mode" && i + 1 < argc) {
            opt.mode = argv[++i];
        } else if (arg == "--batch") {
            opt.batch = true;
        } else if (arg == "--output" && i + 1 < argc) {
            opt.output_dir = argv[++i];
        } else {
            positional.push_back(arg);
        }
    }

    if (opt.mode.empty()) return false;
    if (opt.batch) {
        if (positional.size() < 2) return false;
        opt.vis_dir = positional[0];
        opt.ir_dir = positional[1];
        if (opt.output_dir.empty()) opt.output_dir = "build/output";
    } else {
        if (positional.size() < 2) return false;
        opt.vis_path = positional[0];
        opt.ir_path = positional[1];
        if (opt.output_dir.empty()) opt.output_dir = "build/output";
    }
    return true;
}

int main(int argc, char** argv)
{
    // 使用方式：
    // 單張: FinalProject.exe --mode dwt <vis_path> <ir_path> --output <out_dir> 
    // 批次: FinalProject.exe --mode dwt --batch <vis_dir> <ir_dir> --output <output_dir> 
    Options opt;
    if (!parse_args(argc, argv, opt)) {
        std::cerr << "Wrong arguments.\n"
                  << "  Single: " << argv[0] << " --mode dwt vis.png ir.png --output build/output\n"
                  << "  Batch: " << argv[0] << " --mode dwt --batch assets/Vis assets/Ir --output build/output\n";
        return 1;
    }

    if (opt.mode == "dwt") {
        if (opt.batch) {
            return run_dwt_batch(opt);
        } else {
            return run_dwt_single(opt);
        }
    }

    std::cerr << "[Error] Unsupported mode: " << opt.mode << "\n";
    return 1;
}
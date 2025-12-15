#pragma once

#include <string>

struct Options {
    std::string mode;
    bool batch = false;
    std::string vis_dir;
    std::string ir_dir;
    std::string output_dir = "build/output";
    std::string vis_path;
    std::string ir_path;
};

int run_dwt_single(const Options& opt);
int run_dwt_batch(const Options& opt);

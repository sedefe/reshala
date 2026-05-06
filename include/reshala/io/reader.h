#pragma once

#include "reshala/io/lp/lp_reader.h"
#include "reshala/io/mps/mps_reader.h"

namespace reshala {

reshala::FileReadStatus Read(const char* path, MilpModel& model) {
    std::filesystem::path file_path(path);
    std::string extension = to_lowercase(file_path.extension().string());
    
    if (extension == ".mps") {
        MpsReader mps_reader(file_path, model);
        return mps_reader.Read();
    } else if (extension == ".lp") {
        LpReader lp_reader(file_path, model);
        return lp_reader.Read();
    } else {
        throw std::invalid_argument("Unsupported file format: " + extension);
    }
}

}  // namespace reshala

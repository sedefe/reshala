#pragma once

#include "reshala/io/lp/lp_reader.h"
#include "reshala/io/mps/mps_reader.h"

namespace reshala {

class Io {
   public:
    Io() : mps_reader_(model_, names_), lp_reader_(model_, names_) {}

    reshala::FileReadStatus Read(const char* path) {
        std::filesystem::path file_path(path);
        std::string extension = to_lowercase(file_path.extension().string());

        if (extension == ".mps") {
            return mps_reader_.Read(file_path);
        } else if (extension == ".lp") {
            return lp_reader_.Read(file_path);
        } else {
            throw std::invalid_argument("Unsupported file format: " + extension);
        }
    }

    MilpModel& GetModel() { return model_; }

   private:
    MilpModel model_;
    Names names_;
    MpsReader mps_reader_;
    LpReader lp_reader_;
};

}  // namespace reshala

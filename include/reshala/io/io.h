#pragma once

#include "reshala/io/lp/lp_reader.h"
#include "reshala/io/mps/mps_reader.h"
#include "reshala/reshala.h"

namespace reshala {

class Io {
   public:
    Io(Reshala& ctx) : ctx_(ctx), model_(ctx.GetModel()) {}

    reshala::FileReadStatus Read(const char* path) {
        std::filesystem::path file_path(path);
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Can't open file " << std::filesystem::absolute(path) << "\n";
            return FileReadStatus::kFsError;
        }

        std::string extension = to_lowercase(file_path.extension().string());

        if (extension == ".mps") {
            reader_ = std::make_unique<MpsReader>(model_, names_);
        } else if (extension == ".lp") {
            reader_ = std::make_unique<LpReader>(model_, names_);
        } else {
            std::cerr << "Unsupported file format: " << extension + "\n";
            return FileReadStatus::kFsError;
        }
        try {
            return reader_->Parse(file);
        } catch (const std::runtime_error& e) {
            std::cerr << e.what();
            return FileReadStatus::kParseError;
        }
    }

    const Names& GetNames() const { return names_; }

   private:
    Reshala& ctx_;

    MilpModel& model_;
    Names names_;

    std::unique_ptr<AbstractReader> reader_;
};

}  // namespace reshala

#pragma once

#include "reshala/io/lp/lp_reader.h"
#include "reshala/io/mps/mps_reader.h"

namespace reshala {

class Io {
   public:
    Io() : mps_reader_(model_, names_), lp_reader_(model_, names_) {}

    reshala::FileReadStatus Read(const char* path) {
        std::filesystem::path file_path(path);
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "Can't open file " << std::filesystem::absolute(path) << "\n";
            return FileReadStatus::kFsError;
        }

        std::string extension = to_lowercase(file_path.extension().string());

        if (extension == ".mps") {
            try {
                return mps_reader_.Parse(file);
            } catch (const std::runtime_error& e) {
                std::cerr << e.what();
                return FileReadStatus::kParseError;
            }
        } else if (extension == ".lp") {
            try {
                return lp_reader_.Parse(file);
            } catch (const std::runtime_error& e) {
                std::cerr << e.what();
                return FileReadStatus::kParseError;
            }
        } else {
            std::cerr << "Unsupported file format: " << extension + "\n";
            return FileReadStatus::kFsError;
        }
    }

    MilpModel& GetModel() { return model_; }

    void PrintValues(std::ostream& os, const std::vector<Scalar>& x) const {
        assert(names_.vars.Size() == x.size());
        for (Index iv = 0; iv < x.size(); iv++) {
            if (!IsZero(x[iv])) {
                os << names_.vars.index_to_name[iv] << ": " << x[iv] << " ";
            }
        }
        os << "\n";
    }

   private:
    MilpModel model_;
    Names names_;
    MpsReader mps_reader_;
    LpReader lp_reader_;
};

}  // namespace reshala

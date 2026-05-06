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
            throw std::invalid_argument("Unsupported file format: \"" + extension + "\"");
        }
    }

    MilpModel& GetModel() { return model_; }

    void PrintValues(std::ostream &os, const std::vector<Scalar>& x) const {
        assert(model_.GetNVars() == x.size());
        for (Index iv=0; iv<x.size(); iv++) {
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

#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "reshala/io/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

enum class MpsParseState { kNam, kRow, kCol, kRhs, kRng, kBnd, kDon, kNon };
ExpType MpsChar2ExpType(char c) {
    switch (c) {
        case 'L':
            return ExpType::kLe;
        case 'G':
            return ExpType::kGe;
        case 'E':
            return ExpType::kEq;
        case 'N':
            return ExpType::kNon;
        default:
            assert(false);
            return ExpType::kNon;
    }
}

class MpsReader {
   public:
    MpsReader() {}

    FileReadStatus Read(const char* fname);
    MilpModel& GetModel() { return model; }

   private:
    bool int_marker = false;
    MilpModel model;
    NameMapper var_names;
    NameMapper con_names;
    std::string obj_name;

    std::vector<ExpType> con_types;
    std::vector<Scalar> con_rhs;

    void ParseRows(const std::vector<std::string>&);
    void ParseColumns(const std::vector<std::string>&);
    void ParseRhs(const std::vector<std::string>&);
    void ParseBounds(const std::vector<std::string>&);

    void FinalizeRhs();
};

}  // namespace reshala

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
            throw std::invalid_argument("Unsupported MPS expression type: " + c);
            return ExpType::kNon;
    }
}
enum class MpsBoundType { kLO, kUP, kFX, kFR, kMI, kPL, kBV, kLI, kUI, kNon };
MpsBoundType Str2MpsBoundType(const std::string&);


class MpsReader {
   public:
    MpsReader(const std::filesystem::path& path, MilpModel& model) : path_(path), model_(model) {}
    FileReadStatus Read();

   private:
    const std::filesystem::path& path_;
    MilpModel& model_;
    Index line_number = 0;

    NameMapper var_names;
    NameMapper con_names;
    std::string obj_name;
    bool int_marker = false;
    std::vector<ExpType> con_types;
    std::vector<Scalar> con_rhs;

    void ThrowParseError(const std::string& message) {
        throw std::runtime_error("Line " + std::to_string(line_number) + ": " + message);
    }

    void ParseRows(const std::vector<std::string>&);
    void ParseColumns(const std::vector<std::string>&);
    void ParseRhs(const std::vector<std::string>&);
    void ParseBounds(const std::vector<std::string>&);

    void FinalizeRhs();
};

}  // namespace reshala

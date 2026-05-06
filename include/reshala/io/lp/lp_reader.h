#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "reshala/io/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

enum class LpParseState { kObj, kCon, kBnd, kBin, kGen, kDon, kNon };
ExpType LpChar2ExpType(char c) {
    switch (c) {
        case '<':
            return ExpType::kLe;
        case '>':
            return ExpType::kGe;
        case '=':
            return ExpType::kEq;
        default:
            throw std::invalid_argument("Unsupported LP expression type: " + c);
            return ExpType::kNon;
    }
}

struct Monom {
    Scalar coeff;
    size_t index;
};

class LpReader {
   public:
    LpReader(const std::filesystem::path& path, MilpModel& model) : path_(path), model_(model) {}
    FileReadStatus Read();

   private:
    const std::filesystem::path& path_;
    MilpModel& model_;
    Index line_number = 0;

    NameMapper var_names;

    void ThrowParseError(const std::string& message) {
        throw std::runtime_error("Line " + std::to_string(line_number) + ": " + message);
    }

    void ParseObjective(const std::vector<std::string>&);
    void ParseConstraint(const std::vector<std::string>&);
    void ParseBounds(const std::vector<std::string>&);
    void ParseBinaries(const std::vector<std::string>&);
    void ParseGenerals(const std::vector<std::string>&);

    void ParseLincomb(const std::vector<std::string>& tokens, std::vector<Monom>& lhs,
                      size_t n = -1);
};

}  // namespace reshala

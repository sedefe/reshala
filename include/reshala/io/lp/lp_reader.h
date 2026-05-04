#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "reshala/io/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

enum class LpReadResult { kOk, kFsError, kParseError };
enum class ParseState { kObj, kCon, kBnd, kBin, kGen, kDon, kNon };
enum class ExpType { kGe, kLe, kEq, kNon };
ExpType char2exptype(char c) {
    switch (c) {
        case '<':
            return ExpType::kLe;
        case '>':
            return ExpType::kGe;
        case '=':
            return ExpType::kEq;
        default:
            assert(false);
            return ExpType::kNon;
    }
}

struct Monom {
    Scalar coeff;
    size_t index;
};

class LpReader {
   public:
    LpReader() {}

    LpReadResult Read(const char* fname);
    MilpModel& GetModel() { return model; }

   private:
    MilpModel model;
    NameMapper var_names;

    void ParseObjective(const std::vector<std::string>&);
    void ParseConstraint(const std::vector<std::string>&);
    void ParseBounds(const std::vector<std::string>&);
    void ParseBinaries(const std::vector<std::string>&);
    void ParseGenerals(const std::vector<std::string>&);

    void ParseLincomb(const std::vector<std::string>& tokens, std::vector<Monom>& lhs,
                       size_t n = -1);
};

}  // namespace reshala

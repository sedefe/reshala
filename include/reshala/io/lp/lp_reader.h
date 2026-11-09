#pragma once

#include "reshala/io/abstract_reader.h"
#include "reshala/io/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

enum class LpParseState { kObj, kCon, kBnd, kBin, kGen, kDon, kNon };
inline ExpType LpChar2ExpType(char c) {
    switch (c) {
        case '<':
            return ExpType::kLe;
        case '>':
            return ExpType::kGe;
        case '=':
            return ExpType::kEq;
        default:
            throw std::invalid_argument("Unsupported LP expression type: " + std::string(1, c));
            return ExpType::kNon;
    }
}

struct Monom {
    Scalar coeff;
    Index index;
};

class LpReader : public AbstractReader {
   public:
    LpReader(MilpModel& model, Names& names) : AbstractReader(model, names) {}
    FileReadStatus Parse(std::ifstream& file);

   private:
    bool matrix_finalized = false;
    std::vector<std::string> multiline;

    void ParseObjective(const std::vector<std::string>&);
    void ParseConstraint(const std::vector<std::string>&);
    void ParseBounds(const std::vector<std::string>&);
    void ParseBinaries(const std::vector<std::string>&);
    void ParseGenerals(const std::vector<std::string>&);

    void ParseLincomb(const std::vector<std::string>& tokens, std::vector<Monom>& lhs,
                      Scalar& free_term, Index begin, Index end);
    void FinalizeMatrix();
};

}  // namespace reshala

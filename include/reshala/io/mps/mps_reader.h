#pragma once

#include <unordered_set>

#include "reshala/io/abstract_reader.h"
#include "reshala/io/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

enum class MpsParseState { kNam, kRow, kCol, kRhs, kRng, kBnd, kDon, kNon };
inline ExpType MpsChar2ExpType(char c) {
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
            throw std::invalid_argument("Unsupported MPS expression type: " + std::string(1, c));
            return ExpType::kNon;
    }
}
enum class MpsBoundType { kLO, kUP, kFX, kFR, kMI, kPL, kBV, kLI, kUI, kNon };
MpsBoundType Str2MpsBoundType(const std::string&);

class MpsReader : public AbstractReader {
   public:
    MpsReader(MilpModel& model, Names& names) : AbstractReader(model, names) {}
    FileReadStatus Parse(std::ifstream& file);

   private:
    bool has_obj = false;
    bool int_marker = false;

    std::vector<ExpType> con_types;
    std::unordered_set<std::string> discarded_free_rows;

    void ParseRows(const std::vector<std::string>&);
    void ParseColumns(const std::vector<std::string>&);
    void ParseRhs(const std::vector<std::string>&);
    void ParseRanges(const std::vector<std::string>&);
    void ParseBounds(const std::vector<std::string>&);
};

}  // namespace reshala

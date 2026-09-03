#pragma once

namespace reshala {

enum class RuleType { kTrivial, kFast, kMedium, kExhaustive, kUnknown };
inline RuleType NextLevel(RuleType type, RuleType max_level) {
    if (type == max_level) return RuleType::kUnknown;

    switch (type) {
        case RuleType::kFast:
            return RuleType::kMedium;
        case RuleType::kMedium:
            return RuleType::kExhaustive;
        default:
            return RuleType::kUnknown;
    }
}

enum class RuleResult { kSkipped, kUnchanged, kReduced, kInfeasible, kUnknown };

struct PresolveStat {
    Index n_rm_con = 0;
    Index n_rm_var = 0;
    Index n_ch_bnd = 0;
    Index n_ch_rhs = 0;
    Index n_ch_coeff = 0;

    PresolveStat operator-(const PresolveStat& other) const {
        PresolveStat result;
        result.n_rm_con = n_rm_con - other.n_rm_con;
        result.n_rm_var = n_rm_var - other.n_rm_var;
        result.n_ch_bnd = n_ch_bnd - other.n_ch_bnd;
        result.n_ch_rhs = n_ch_rhs - other.n_ch_rhs;
        result.n_ch_coeff = n_ch_coeff - other.n_ch_coeff;
        return result;
    }
};

class BitMask {
    std::vector<uint64_t> data;

   public:
    explicit BitMask(Index size) : data((size + 63) / 64, 0) {}
    void Set(Index pos) { data[pos / 64] |= (1ULL << (pos % 64)); }
    bool Get(Index pos) const { return (data[pos / 64] >> (pos % 64)) & 1; }

    void Clear() { std::fill(data.begin(), data.end(), 0); }
};

}  // namespace reshala

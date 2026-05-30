#pragma once

namespace reshala {

enum class RuleType { kTrivial, kFast, kMedium, kExhaustive, kUnknown };
inline RuleType NextLevel(RuleType type) {
    switch (type) {
        case RuleType::kFast:
            return RuleType::kMedium;
        case RuleType::kMedium:
            return RuleType::kExhaustive;
        default:
            return RuleType::kUnknown;
    }
}

enum class RuleResult { kSkipped, kUnchanged, kReduced, kUnknown };

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

}  // namespace reshala

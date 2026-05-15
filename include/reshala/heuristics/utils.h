#pragma once

#include "reshala/model/milp_model.h"

namespace reshala {

enum class RoundingType {
    kAll = 0,
    kInts = 1,
    kNone = 2,
};
std::string RoundingType2Str(RoundingType type);

void Fixing(RoundingType type, MilpModel &model, const std::vector<Scalar> &relaxed_x);

}  // namespace reshala

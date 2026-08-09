#pragma once

#include "reshala/model/milp_model.h"

namespace reshala {

enum class FixingType {
    kAll = 0,
    kInts = 1,
    kNone = 2,
};
std::string FixingType2Str(FixingType type);

void Fixing(FixingType type, MilpModel &model, const std::vector<Scalar> &relaxed_x);

}  // namespace reshala

#pragma once

#include "reshala/presolve/activity.h"
#include "reshala/presolve/utils.h"

namespace reshala {

class ModelTracker;

struct Bounder {
    Bounder(const ModelTracker& t);

    void Reset();

    const ModelTracker& tracker;

    std::vector<Activity> activities;
    Domain domain;

    std::vector<Index> changed_cons;
    BitMask con_mask;
    std::vector<Index> changed_vars;
    BitMask var_mask;

    bool Propagate(Index iv_start, const Bounds& new_bnd);
};

}  // namespace reshala

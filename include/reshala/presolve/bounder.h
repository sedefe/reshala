#pragma once

#include "reshala/presolve/activity.h"
#include "reshala/presolve/utils.h"

namespace reshala {

class ModelTracker;

struct Bounder {
    Bounder(const ModelTracker& t);

    void Reset();
    bool Propagate(Index iv_start, const Bounds& new_bnd);

    const ModelTracker& tracker;

    std::vector<Activity> activities;
    Domain domain;

    std::vector<Index> changed_cons;
    BitMask changed_con_mask;
    std::vector<Index> changed_vars;
    BitMask changed_var_mask;

    std::vector<Index> all_changed_cons;
    BitMask all_changed_con_mask;
    std::vector<Index> all_changed_vars;
    BitMask all_changed_var_mask;

    BitMask redundant_con_mask;
};

}  // namespace reshala

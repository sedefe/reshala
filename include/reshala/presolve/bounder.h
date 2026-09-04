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

    MaskedVector changed_cons;
    MaskedVector changed_vars;

    MaskedVector all_changed_vars;
    MaskedVector redundant_con_mask;
};

}  // namespace reshala

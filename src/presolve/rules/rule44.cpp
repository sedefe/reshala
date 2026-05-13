#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule44::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    // Todo: count locks?
    // Обратный порядок тут нужен, чтобы в постсолве индексация не поехала
    for (Index iv = model.GetNVars() - 1; iv >= 0; iv--) {
        if (info.GetVarMask(iv)) continue;

        bool eligible_up = model.GetObj().coefficients[iv] <= 0;
        bool eligible_down = model.GetObj().coefficients[iv] >= 0;

        for (SvIterator el(model.GetAc().GetCol(iv)); el & (eligible_up || eligible_down); ++el) {
            const Bounds& rhs = model.GetRhs()[el.index()];
            if (rhs.le == -kInf) {
                eligible_down &= el.value() >= 0;
                eligible_up &= el.value() <= 0;
            }
            if (rhs.ri == kInf) {
                eligible_down &= el.value() <= 0;
                eligible_up &= el.value() >= 0;
            }
            if (rhs.le != -kInf or rhs.ri != kInf) {
                eligible_down = eligible_up = false;
            }
        }

        if (eligible_down) {
            FixVariableTransform tr(iv, model.GetBounds(iv).le);
            tr.Do(info);
            transforms.push_back(std::make_unique<FixVariableTransform>(tr));

            n_reduced++;
        } else if (eligible_up) {
            FixVariableTransform tr(iv, model.GetBounds(iv).ri);
            tr.Do(info);
            transforms.push_back(std::make_unique<FixVariableTransform>(tr));

            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala

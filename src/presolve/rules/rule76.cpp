#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule76::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    // Primal
    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const Bounds& rhs = model.GetRhs(ic);
        if (StrongLt(rhs.le, rhs.ri)) continue;
        Scalar b = (rhs.le + rhs.ri) / 2;
        if (!IsZero(MinFraction(b))) continue;

        const auto& row = model.GetRow(ic);

        Index n_continuous = 0;
        Index iv_cont = -1;
        Scalar a = kNan;
        bool eligible = true;
        for (SvIterator el(row); el and eligible; ++el) {
            if (tracker.GetVarMask(el.index())) continue;

            if (model.GetIntegrality(el.index())) {
                eligible &= IsZero(MinFraction(el.value()));
            } else {
                n_continuous++;
                iv_cont = el.index();
                a = el.value();
                eligible &= (n_continuous <= 1);
                eligible &= IsZero(MinFraction(a));
            }
        }
        eligible &= (n_continuous == 1);
        if (eligible) {
            tracker.ScaleVar(iv_cont, 1 / std::abs(a));
            tracker.SetVarInt(iv_cont);
            n_reduced++;
        }
    }

    // Todo: Dual

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala

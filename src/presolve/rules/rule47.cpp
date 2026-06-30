#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule47::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;
    Index size = 0;

    // {ic: [(iv1, a1), (iv2, a2), ...]}
    std::unordered_map<Index, std::vector<std::pair<Index, Scalar>>> candidates;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.GetVarMask(iv)) continue;

        Index ic = -1;
        Scalar a = kNan;
        for (SvIterator el(model.GetCol(iv)); el; ++el) {
            if (tracker.GetConMask(el.index())) continue;
            if (ic != -1) {
                ic = -1;
                break;
            }
            ic = el.index();
            a = el.value();
        }
        if (ic < 0) continue;
        const auto& rhs = model.GetRhs(ic);
        if (!IsZero(rhs.ri - rhs.le)) continue;
        if (model.GetIntegrality(iv)) {  // Handle integer-row slacks with |a|=1
            if (!IsZero(std::abs(a) - 1)) continue;
            if (!model.RowIsInteger(ic)) continue;
            candidates[ic].push_back({iv, a});
        } else {
            candidates[ic].push_back({iv, a});
        }
    }

    // Если слаков несколько, выбираем тот, у которого меньше c[iv]
    // Todo Должна быть нормальная логика
    for (const auto& [ic, vec] : candidates) {
        Index iv_cand;
        Scalar a_cand;
        Scalar min_score = kInf;
        for (const auto& [iv, a] : vec) {
            Scalar score = model.GetObj().coefficients[iv];
            if (score < min_score) {
                min_score = score;
                iv_cand = iv;
                a_cand = a;
            }
        }
        tracker.SlackSub(ic, iv_cand, a_cand);
        n_reduced++;
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala

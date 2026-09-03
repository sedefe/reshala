#include <array>

#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule72::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;
    tracker.GetImplications().clear();

    std::array<Bounder, 2> bounders = {Bounder{tracker}, Bounder{tracker}};
    std::array<bool, 2> results;

    std::vector<std::pair<Index, Index>> scored_binaries;
    scored_binaries.reserve(model.GetNVars());
    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.GetVarMask(iv)) continue;
        if (!model.IsBinary(iv)) continue;
        Index score = 0;
        for (SvIterator el(model.GetCol(iv)); el; ++el) {
            score += model.GetRow(el.index()).Size() == 2;  // likely an implication
        }
        scored_binaries.push_back({-score, iv});
    }
    std::sort(scored_binaries.begin(), scored_binaries.end());

    for (const auto& pair : scored_binaries) {
        Index iv = pair.second;
        if (tracker.GetVarMask(iv)) continue;
        if (!model.IsBinary(iv)) continue;

        for (Index i = 0; i < 2; i++) {
            bounders[i].Reset();
            results[i] = bounders[i].Propagate(iv, {Scalar(i), Scalar(i)});
        }

        if (!results[0] and !results[1]) {
            return RuleResult::kInfeasible;
        }
        if (!results[0]) {
            tracker.ImportBounder(bounders[1]);
            n_reduced++;
            continue;
        }
        if (!results[1]) {
            tracker.ImportBounder(bounders[0]);
            n_reduced++;
            continue;
        }

        for (Index iv1 = 0; iv1 < model.GetNVars(); iv1++) {
            if (tracker.GetVarMask(iv1)) continue;
            if (iv == iv1) continue;
            const Bounds& bnd = model.GetBounds(iv1);
            const Bounds& bnd0 = bounders[0].domain.GetBounds(iv1);
            const Bounds& bnd1 = bounders[1].domain.GetBounds(iv1);
            Bounds derived = {
                std::min(bnd0.le, bnd1.le),
                std::max(bnd0.ri, bnd1.ri),
            };

            // Check if we can fix or substitute x <- a*y + b
            if (WeakEq(bnd0.le, bnd0.ri) and WeakEq(bnd1.le, bnd1.ri)) {
                Scalar y0 = (bnd0.le + bnd0.ri) / 2;
                Scalar y1 = (bnd1.le + bnd1.ri) / 2;
                if (WeakEq(y0, y1)) {
                    tracker.FixVar(iv1, (y0 + y1) / 2);
                } else {
                    tracker.SimpleSub(iv1, y1 - y0, iv, y0);
                }
                n_reduced++;
                continue;
            }

            // Check if we can strengthen the bounds
            if (StrongGt(derived.le, bnd.le) or StrongLt(derived.ri, bnd.ri)) {
                Bounds new_bnd = {std::max(bnd.le, derived.le), std::min(bnd.ri, derived.ri)};
                if (StrongGt(new_bnd.le, new_bnd.ri)) {
                    return RuleResult::kInfeasible;
                }
                tracker.UpdVarBounds(iv1, std::move(new_bnd));
                n_reduced++;
                continue;
            }

            if (StrongGt(bnd0.le, bnd.le)) {  // x = 0 => y >= b
                tracker.GetImplications().push_back({iv, false, iv1, true, bnd0.le});
            }
            if (StrongLt(bnd0.ri, bnd.ri)) {  // x = 0 => y <= b
                tracker.GetImplications().push_back({iv, false, iv1, false, bnd0.ri});
            }
            if (StrongGt(bnd1.le, bnd.le)) {  // x = 1 => y >= b
                tracker.GetImplications().push_back({iv, true, iv1, true, bnd1.le});
            }
            if (StrongLt(bnd1.ri, bnd.ri)) {  // x = 1 => y <= b
                tracker.GetImplications().push_back({iv, true, iv1, false, bnd1.ri});
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala

#include <array>

#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule72::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    tracker.GetImplications().clear();

    Index n_reduced = 0, n_probed = 0;
    const Index kMaxReductions = std::max(1000, (model.GetNCons() + model.GetNVars()) / 20);

    std::array<Bounder, 2> bounders = {Bounder{tracker}, Bounder{tracker}};
    std::array<bool, 2> results;

    Index find_cliq_res = FindCliques(tracker);
    n_reduced += find_cliq_res;
    Score(tracker);

    for (const auto& cand : scores) {
        Index iv = cand[3];
        if (tracker.VarIsDeleted(iv)) continue;
        if (!model.IsBinary(iv)) continue;

        n_probed++;
        if (n_reduced >= kMaxReductions) break;

        for (Index i = 0; i < 2; i++) {
            bounders[i].Reset();
            results[i] = bounders[i].Propagate(iv, {Scalar(i), Scalar(i)});
        }

        if (!results[0] and !results[1]) {
            return RuleResult::kInfeasible;
        }
        if (!results[0]) {
            n_reduced += tracker.ImportBounder(bounders[1]);
            continue;
        }
        if (!results[1]) {
            n_reduced += tracker.ImportBounder(bounders[0]);
            continue;
        }

        for (Index iv1 = 0; iv1 < model.GetNVars(); iv1++) {
            if (tracker.VarIsDeleted(iv1)) continue;
            if (iv == iv1) continue;
            const Bounds& bnd = model.GetBounds(iv1);
            const Bounds& bnd0 = bounders[0].domain.GetBounds(iv1);
            const Bounds& bnd1 = bounders[1].domain.GetBounds(iv1);

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
            Bounds derived = {
                std::min(bnd0.le, bnd1.le),
                std::max(bnd0.ri, bnd1.ri),
            };
            assert(StrongLt(derived.le, derived.ri));
            if (GoodBndChange(bnd, derived)) {
                tracker.UpdVarBounds(iv1, std::move(derived));
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

Index Rule72::FindCliques(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();

    Index res = 0;
    cliq_implics.clear();

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        const auto& row = model.GetRow(ic);
        const auto& rhs = model.GetRhs(ic);

        auto try_clique = [&](const SparseVector& row, Scalar b, bool negate) {
            if (b == kInf) return;
            Index n_neg = 0;
            for (SvIterator el(row); el; ++el) {
                if (tracker.VarIsDeleted(el.index())) continue;
                if (!model.IsBinary(el.index())) {
                    return;
                }
                Scalar a = negate ? -el.value() : el.value();
                if (a == 1.0) {
                } else if (a == -1.0)
                    ++n_neg;
                else
                    return;
            }
            if (StrongLt(b + n_neg, 1.0)) return;
            if (WeakGe(b + n_neg, 2.0)) return;

            for (SvIterator el(row); el; ++el) {
                bool val = (el.value() > 0) != negate;
                auto [iter, _] = cliq_implics.try_emplace({el.index(), val});
                iter->second.reserve(iter->second.size() + row.Size() - 1);  // <len(row) emplacings
                for (SvIterator el1(row); el1; ++el1) {
                    if (el.index() == el1.index()) continue;
                    bool val1 = (el1.value() > 0) != negate;
                    iter->second.emplace(BinValue{el1.index(), !val1});
                }
            }
        };

        // if (row.Size() > 2) continue;
        try_clique(row, rhs.ri, false);  // <= b
        try_clique(row, -rhs.le, true);  // >= b, negated
    }

    if (0) {  // Если в импликациях есть пересечения, можно форсить
        for (Index iv = 0; iv < model.GetNVars(); iv++) {
            if (tracker.VarIsDeleted(iv)) continue;
            if (!model.IsBinary(iv)) continue;

            auto it0 = cliq_implics.find({iv, 0});
            auto it1 = cliq_implics.find({iv, 1});
            if (it0 != cliq_implics.end() and it1 != cliq_implics.end()) {
                for (const auto& v : it0->second) {
                    if (tracker.VarIsDeleted(v.iv)) continue;
                    auto it2 = it1->second.find(v);
                    if (it2 != it1->second.end()) {
                        tracker.FixVar(v.iv, v.val);
                        res++;
                    }
                }
            }
        }
    }

    return res;
}

void Rule72::Score(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();

    scores.clear();
    scores.reserve(model.GetNVars());
    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.VarIsDeleted(iv)) continue;
        if (!model.IsBinary(iv)) continue;
        Index impl_d = cliq_implics[{iv, 0}].size();
        Index impl_u = cliq_implics[{iv, 1}].size();

        Index score1 = std::min(5000, impl_d * impl_u);
        Index score2 = std::min(100, impl_d + impl_u);
        Index score3 = rand();

        scores.push_back({-score1, -score2, score3, iv});
    }
    std::sort(scores.begin(), scores.end());
}

}  // namespace reshala

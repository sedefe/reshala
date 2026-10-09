#include <unordered_map>

#include "reshala/presolve/rules.h"

namespace reshala {

void SortOffsets(std::vector<Index>& target, const std::vector<Index>& indices, const Locks& locks,
                 size_t max_len) {
    if (max_len <= 0 || target.empty()) return;

    std::size_t k = std::min(target.size(), max_len);
    std::partial_sort(target.begin(), target.begin() + k, target.end(), [&](Index a, Index b) {
        return locks.sum_locks[indices[a]] < locks.sum_locks[indices[b]];
    });
}

RuleResult Rule53::Apply(ModelTracker& tracker) {
    Index n_reduced = 0;

    tracker.InitLocks();  // Todo: keep & update

    FillEqHashMap(tracker);
    if (eq_hash_map.size() > 0) {
        n_reduced = PairSearch(tracker);
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

void Rule53::FillEqHashMap(const ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();

    eq_hash_map.clear();
    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const Bounds& rhs = model.GetRhs(ic);
        if (!WeakEq(rhs.le, rhs.ri)) continue;

        const auto& row = model.GetRow(ic);

        std::vector<Index> indices(row.Size());
        std::iota(indices.begin(), indices.end(), 0);
        SortOffsets(indices, row.indices(), model.GetLocks(), kMaxNzs);
        Index sz = std::min(row.Size(), kMaxNzs);

        for (Index nz2 = 0; nz2 < sz; ++nz2) {
            Index iv2 = row.indices()[indices[nz2]];
            Scalar a2 = row.values()[indices[nz2]];
            for (Index nz1 = 0; nz1 < nz2; ++nz1) {
                Index iv1 = row.indices()[indices[nz1]];
                Scalar a1 = row.values()[indices[nz1]];
                const Key key{iv1, iv2, a2 / a1};
                const Value value{ic, a1, a2};

                auto [it, inserted] = eq_hash_map.try_emplace(key, value);
                Value& v = it->second;
                if (!inserted) {
                    if (row.Size() < model.GetRow(v.ic).Size()) v = value;
                }
            }
        }
    }
}

Index Rule53::PairSearch(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_paired = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const auto& row = model.GetRow(ic);

        while (true) {
            Index min_size = row.Size();
            Index ic_pair = -1;
            Scalar best_lambda;

            std::vector<Index> offsets(row.Size());
            std::iota(offsets.begin(), offsets.end(), 0);
            SortOffsets(offsets, row.indices(), model.GetLocks(), kMaxNzs);
            Index sz = std::min(row.Size(), kMaxNzs);

            for (Index nz2 = 0; nz2 < sz; ++nz2) {
                Index iv2 = row.indices()[offsets[nz2]];
                Scalar a2 = row.values()[offsets[nz2]];
                for (Index nz1 = 0; nz1 < nz2; ++nz1) {
                    Index iv1 = row.indices()[offsets[nz1]];
                    Scalar a1 = row.values()[offsets[nz1]];
                    const Key key{iv1, iv2, a2 / a1};

                    auto it = eq_hash_map.find(key);
                    if (it != eq_hash_map.end()) {
                        const Value& value = it->second;
                        if (value.ic == ic) continue;

                        Scalar lambda = a1 / value.a1;
                        Index size = axpy_size(-lambda, model.GetRow(value.ic), row);

                        if (size < min_size) {
                            min_size = size;
                            ic_pair = value.ic;
                            best_lambda = lambda;
                            if (min_size == 0) break;
                        }
                    }
                }
            }

            if (ic_pair >= 0) {
                n_paired++;

                SparseVector new_row = axpy(-best_lambda, model.GetRow(ic_pair), row);

                const Bounds& rhs = model.GetRhs(ic);
                Scalar b = (model.GetRhs(ic_pair).le + model.GetRhs(ic_pair).ri) / 2;
                tracker.UpdCon(ic, new_row);
                tracker.UpdActivity(ic);
                tracker.UpdRhs(ic, {rhs.le - best_lambda * b, rhs.ri - best_lambda * b});
            } else {
                break;
            }
        }
    }

    if (n_paired > 0) {
        tracker.RebuildAc();
    }

    return n_paired;
}

}  // namespace reshala

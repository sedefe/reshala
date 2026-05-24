#include <unordered_map>

#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule52::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    std::unordered_map<Index, std::vector<Index>> bins;
    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (info.GetConMask(ic)) continue;
        if (model.GetAr().GetRow(ic).Empty()) continue;
        auto hash = HashRow(info, ic);

        bins[hash].push_back(ic);
    }

    for (const auto& pair : bins) {
        const auto& bin = pair.second;
        if (bin.size() == 1) continue;
        std::vector<std::vector<Index>> sub_bins;

        bin_scales.assign(bin.size(), 0);
        for (Index i = 0; i < bin.size(); i++) {
            for (const auto& value : model.GetAr().GetRow(bin[i]).values()) {
                if (StrongGt(std::abs(value), std::abs(bin_scales[i]))) {
                    // Тут StrongGt, чтобы среди одинаково больших коэффициентов выбрать первый.
                    // Если строки параллельные, то максимальные элементы окажутся на одних местах.
                    bin_scales[i] = value;
                }
            }
        }

        sub_bins.push_back({0});
        for (Index i = 1; i < bin.size(); i++) {
            bool not_par = true;
            for (auto& sub_bin : sub_bins) {
                if (Parallel(model.GetAr().GetRow(bin[sub_bin[0]]), bin_scales[sub_bin[0]],
                             model.GetAr().GetRow(bin[i]), bin_scales[i])) {
                    sub_bin.push_back(i);
                    not_par = false;
                    break;
                }
            }
            if (not_par) {
                sub_bins.push_back({bin[i]});
            }
        }
        for (auto& sub_bin : sub_bins) {
            if (sub_bin.size() == 1) continue;
            Scalar max_abs_scale = 0;
            Index i_max_abs = -1;
            for (auto i : sub_bin) {
                if (std::abs(bin_scales[i]) > std::abs(max_abs_scale)) {
                    i_max_abs = i;
                    max_abs_scale = bin_scales[i];
                }
            }

            Bounds& rhs = model.GetRhs(bin[i_max_abs]);
            for (Index i : sub_bin) {
                if (i == i_max_abs) continue;
                info.MaskCon(bin[i]);
                Bounds rhs1 = model.GetRhs(bin[i]);
                auto scale = bin_scales[i];
                if (std::signbit(scale) == std::signbit(max_abs_scale)) {
                    rhs1 = {rhs1.le * max_abs_scale / scale, rhs1.ri * max_abs_scale / scale};
                } else {
                    rhs1 = {rhs1.ri * max_abs_scale / scale, rhs1.le * max_abs_scale / scale};
                }

                rhs = BoundsIntersection(rhs, rhs1);
                n_reduced++;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

Index Rule52::HashRow(const ModelInfo& info, Index ic) const {
    const auto& row = info.GetModel().GetAr().GetRow(ic);
    Index hash = row.Size();
    for (Index iv : row.indices()) {
        if (info.GetVarMask(iv)) continue;
        hash = hash * 0x9E3779B9 + (iv + 0x3AB0D);
    }
    return hash;
}

bool Rule52::Parallel(const SparseVector& sv1, Scalar scale1, const SparseVector& sv2,
                      Scalar scale2) const {
    if (sv1.Size() != sv2.Size()) return false;
    if (sv1.Empty()) return true;
    if (sv1.indices() != sv2.indices()) return false;

    Index n = sv1.Size();
    const Scalar* p1 = sv1.values().data();
    const Scalar* p2 = sv2.values().data();
    Scalar mul = scale2 / scale1;
    if (std::abs(scale1) < std::abs(scale2)) {
        std::swap(p1, p2);
        mul = scale1 / scale2;
    }

    for (Index i = 0; i < n; i++) {
        if (!IsZero(p1[i] * mul - p2[i])) return false;
    }

    return true;
}

}  // namespace reshala

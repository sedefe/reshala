#pragma once

#include "reshala/model/milp_model.h"

namespace reshala {

enum class CutType { kProbing = 0, kCmir = 1 };
std::string CutType2Str(CutType type);

struct Cut {
    // lhs >= rhs

    Cut(CutType t, const SparseVector& l, Scalar r) : type(t), lhs(l), rhs(r) {}

    SparseVector lhs;
    Scalar rhs;

    Scalar quality = kNan;

    bool removed = false;
    bool selected = false;
    CutType type;
    Index age = 0;
    Index round = -1;
    Scalar norm2 = kNan;

    Scalar m_violation = kNan;
    Scalar m_rel_violation = kNan;
    Scalar m_distance = kNan;
    Scalar m_distance1 = kNan;
    Scalar m_distance2 = kNan;
    Scalar m_obj_parallelism = kNan;
    Scalar m_support = kNan;
    Scalar m_int_support = kNan;
    Scalar m_exp_improvement = kNan;
    Scalar m_unit_act_distance = kNan;

    void CalcMetrics(const Solution& sol, const MilpModel& model) {
        norm2 = lhs.Norm2();

        m_violation = GetViolation(sol.x);
        m_rel_violation = GetRelViolation(sol.x);
        m_obj_parallelism = GetObjParallelism(model.GetObj().coefficients);
        m_distance = GetDistance(sol.x);
        m_distance1 = GetDistanceI(sol.x, model);
        m_distance2 = GetDistanceII(sol.x);
        m_exp_improvement = GetExpectedImprovement();
        m_support = GetSupport();
        m_int_support = GetIntegralSupport();
        m_unit_act_distance = GetUnitActiveDistance(sol.x, model);
    }

    bool IsViolated(const std::vector<Scalar>& x) const { return StrongGt(GetViolation(x), 0.0); }

   private:
    //////////////////////////////////////////////////
    // Basic rule for all metrics: Bigger is better //
    //////////////////////////////////////////////////

    inline Scalar GetViolation(const std::vector<Scalar>& x) const {
        // v(a, b, x) = b - a * x
        Scalar l;
        dot(lhs, x, l);
        return rhs - l;
    }

    inline Scalar GetRelViolation(const std::vector<Scalar>& x) const {
        // d(a, b, x) = v(a, b, x) / max(|b|, 1)
        Scalar denominator = rhs;
        return m_violation / std::max(std::abs(denominator), kEpsZero);
    }

    inline Scalar GetDistance(const std::vector<Scalar>& x) const {
        // d(a, b, x) = v(a, b, x) / ||a||
        // Считаем квадрат этой метрики
        return (m_violation * m_violation) / norm2;
    }

    inline Scalar GetDistanceI(const std::vector<Scalar>& x, const MilpModel& model) const {
        Scalar denominator(1);
        for (SvIterator el(lhs); el; ++el) {
            auto iv = el.index();
            const auto& bnd = model.GetBounds(iv);
            if (StrongGt(x[iv], bnd.le) and StrongLt(x[iv], bnd.ri)) {
                denominator += el.value() * el.value();
            }
        }
        return (m_violation * m_violation) / denominator;
    }

    inline Scalar GetDistanceII(const std::vector<Scalar>& x) const { return kNan; }

    inline Scalar GetObjParallelism(const std::vector<Scalar>& obj) const {
        // o(a) = ac / (||a|| * ||c||)
        // Считаем квадрат этой метрики
        Scalar scal_mult;
        dot(lhs, obj, scal_mult);
        auto cosine2 = (scal_mult * scal_mult) / (norm2 * SparseVector(obj).Norm2());
        return cosine2;
    }

    inline Scalar GetExpectedImprovement() const {
        // e(a, b, x) = ||c|| * o(a) * d(a, b, x)
        // Считаем квадрат этой метрики
        // Игнорируем ||c||, квадраты остальных множителей посчитаны выше
        return m_obj_parallelism * m_distance;
    }

    inline Scalar GetSupport() const { return Scalar(lhs.Size()); }

    inline Scalar GetIntegralSupport() const { return kNan; }

    inline Scalar GetUnitActiveDistance(const std::vector<Scalar>& x,
                                        const MilpModel& model) const {
        // v(a, b, x) / norm_act / n_act_nz
        // Считаем квадрат этой метрики
        Index n_active_nz = 0;
        Scalar norm_active = 0.;

        for (SvIterator el(lhs); el; ++el) {
            auto iv = el.index();
            const auto& bnd = model.GetBounds(iv);
            if (el.value() > 0) {
                if (StrongLt(x[iv], bnd.ri)) {
                    n_active_nz++;
                    norm_active += el.value() * el.value();
                }
            } else {
                if (StrongGt(x[iv], bnd.le)) {
                    n_active_nz++;
                    norm_active += el.value() * el.value();
                }
            }
        }

        if (n_active_nz == 0) return 0.;
        return (m_violation * m_violation) / norm_active / n_active_nz / n_active_nz;
    }
};

std::ostream& operator<<(std::ostream& os, const Cut& cut);

}  // namespace reshala

#pragma once

#include <memory>
#include <unordered_map>

#include "reshala/model/implication.h"
#include "reshala/model/milp_model.h"
#include "reshala/presolve/activity.h"
#include "reshala/presolve/bounder.h"
#include "reshala/presolve/transforms.h"
#include "reshala/presolve/utils.h"

namespace reshala {

class ModelTracker {
   public:
    ModelTracker(MilpModel& model);
    PresolveStat stat;

    const MilpModel& GetModel() const { return model_; }

    void MaskCon(Index ic) {
        deleted_cons_.Add(ic);
        stat.n_rm_con++;
    }
    void MaskVar(Index iv) {
        deleted_vars_.Add(iv);
        stat.n_rm_var++;
    }
    inline bool GetConMask(Index ic) const { return deleted_cons_.Get(ic); }
    inline bool GetVarMask(Index ic) const { return deleted_vars_.Get(ic); }

    void CompressCons();
    void CompressVars();

    void CalcActivities();
    Activity CalcActivity(Index ic) const;
    inline const std::vector<Activity>& GetActivities() const { return activities_; }
    inline const Activity& GetActivity(Index ic) const { return activities_[ic]; }
    inline const Bounds GetConRange(Index ic) const { return activities_[ic].GetRange(); }

    inline std::vector<Implication>& GetImplications() { return implications_; }
    inline const std::vector<Implication>& GetImplications() const { return implications_; }

    inline Index GetNDeletedCons() const { return deleted_cons_.GetNValues(); }
    inline Index GetNDeletedVars() const { return deleted_vars_.GetNValues(); }
    inline const MaskedVector& GetDeletedCons() const { return deleted_cons_; }
    inline const MaskedVector& GetDeletedVars() const { return deleted_vars_; }

    // Model transformations
    void FixVar(Index iv, Scalar val);
    void ConstShiftVar(Index iv, Scalar val);
    bool SimpleSub(Index iv1, Scalar a, Index iv2, Scalar b);  // iv1 <- a*iv2 + b
    void SlackSub(Index ic, Index iv, Scalar a);

    void UpdRhs(Index ic, Bounds rhs);
    void UpdVarBounds(Index iv, Bounds bnd);
    void UpdCoeff(Index ic, Index iv, Scalar val);
    void ScaleObj(Scalar x);
    void ScaleObjExp(Index e);
    void ScaleRow(Index ic, Scalar x);
    void ImportBounder(Bounder& bounder);

    inline Index GetOrigNVars() const { return orig_n_vars_; }
    inline const std::vector<Index>& GetOrigVarIdx() const { return orig_var_idx_; }

    inline const std::vector<std::unique_ptr<Transform>>& GetTransforms() const {
        return transforms_;
    }

    Bounds DeriveBounds(Index ic, Index iv, Activity act, const Bounds& bnd, Scalar val) const;

   private:
    MilpModel& model_;

    Index orig_n_vars_;
    std::vector<Index> orig_var_idx_;

    MaskedVector deleted_cons_;
    MaskedVector deleted_vars_;

    std::vector<Activity> activities_;

    std::vector<Implication> implications_;

    std::vector<std::unique_ptr<Transform>> transforms_;
};

}  // namespace reshala

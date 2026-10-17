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

    // Masks
    inline void MaskCon(Index ic) {
        deleted_cons_.Add(ic);
        stat.n_rm_con++;
    }
    inline void MaskVar(Index iv) {
        deleted_vars_.Add(iv);
        stat.n_rm_var++;
    }
    inline const MaskedVector& GetDeletedCons() const { return deleted_cons_; }
    inline const MaskedVector& GetDeletedVars() const { return deleted_vars_; }
    inline Index GetNDeletedCons() const { return deleted_cons_.GetNValues(); }
    inline Index GetNDeletedVars() const { return deleted_vars_.GetNValues(); }
    inline bool ConIsDeleted(Index ic) const { return deleted_cons_.Get(ic); }
    inline bool VarIsDeleted(Index ic) const { return deleted_vars_.Get(ic); }

    inline const MaskedVector& GetChangedCons() const { return changed_cons_; }
    inline const MaskedVector& GetChangedVars() const { return changed_vars_; }
    inline void TouchAll() {
        changed_cons_.SetAll();
        changed_vars_.SetAll();
    };
    inline void TouchCon(Index ic) { changed_cons_.Add(ic); }
    inline void TouchVar(Index iv) { changed_vars_.Add(iv); }
    inline void TouchConVars(Index ic) {
        for (auto iv : model_.GetRow(ic).indices()) changed_vars_.Add(iv);
    }
    inline void TouchVarCons(Index iv) {
        for (auto ic : model_.GetCol(iv).indices()) changed_cons_.Add(ic);
    }
    inline void ClearChangedCons() { changed_cons_.Clear(); }
    inline void ClearChangedVars() { changed_vars_.Clear(); }

    void CompressCons();
    void CompressVars();

    // Tracking info
    inline const std::vector<Activity>& GetActivities() const { return activities_; }
    inline const Activity& GetActivity(Index ic) const { return activities_[ic]; }
    inline const Bounds GetConRange(Index ic) const { return activities_[ic].GetRange(); }

    inline const Locks& GetLocks() const { return locks_; }
    void CalcLocks();
    void UpdLocks(Index iv);
    void CalcActivities();
    void UpdActivity(Index ic);

    inline std::vector<Implication>& GetImplications() { return implications_; }
    inline const std::vector<Implication>& GetImplications() const { return implications_; }

    // Model transformations
    void FixVar(Index iv, Scalar val);
    void ConstShiftVar(Index iv, Scalar val);
    bool SimpleSub(Index iv1, Scalar a, Index iv2, Scalar b);
    void SlackSub(Index ic, Index iv, Scalar a);
    void UpdCon(Index ic, const SparseVector& row);
    void RebuildAc();

    void UpdRhs(Index ic, Bounds rhs);
    void UpdVarBounds(Index iv, Bounds bnd);
    void UpdCoeff(Index ic, Index iv, Scalar val);
    void ScaleObj(Scalar x);
    void ScaleObjExp(Index e);
    void ScaleCon(Index ic, Scalar scale);
    void ScaleVar(Index iv, Scalar scale);
    void SetVarInt(Index iv);
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

    MaskedVector changed_cons_, changed_vars_;

    MaskedVector deleted_cons_;
    MaskedVector deleted_vars_;

    Locks locks_;
    std::vector<Activity> activities_;

    std::vector<Implication> implications_;

    std::vector<std::unique_ptr<Transform>> transforms_;
};

}  // namespace reshala

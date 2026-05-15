#pragma once

#include "reshala/model/milp_model.h"

namespace reshala {

class BitMask {
    std::vector<uint64_t> data;

   public:
    explicit BitMask(Index size) : data((size + 63) / 64, 0) {}
    void Set(Index pos) { data[pos / 64] |= (1ULL << (pos % 64)); }
    bool Get(Index pos) const { return (data[pos / 64] >> (pos % 64)) & 1; }

    void Clear() { std::fill(data.begin(), data.end(), 0); }
};

class ModelInfo {
   public:
    ModelInfo(MilpModel& model);

    const MilpModel& GetModel() const { return model_; }
    MilpModel& GetModel() { return model_; }

    void MaskCon(Index ic) {
        con_mask_.Set(ic);
        deleted_cons_.push_back(ic);
    }
    void MaskVar(Index iv) {
        var_mask_.Set(iv);
        deleted_vars_.push_back(iv);
    }
    inline bool GetConMask(Index ic) const { return con_mask_.Get(ic); }
    inline bool GetVarMask(Index ic) const { return var_mask_.Get(ic); }

    void CompressCons();
    void CompressVars();

    void CalcActivities();
    const Bounds& GetActivity(Index ic) const { return activities_[ic]; }
    Bounds& GetActivity(Index ic) { return activities_[ic]; }

    Index GetNDeletedCons() const { return deleted_cons_.size(); }
    Index GetNDeletedVars() const { return deleted_vars_.size(); }
    const std::vector<Index>& GetDeletedCons() const { return deleted_cons_; }
    const std::vector<Index>& GetDeletedVars() const { return deleted_vars_; }

    void FixVar(Index iv, Scalar val);
    void UpdVarBounds(Index iv, const Bounds& bnd);

    inline Index GetOrigNVars() const { return orig_n_vars_; }
    inline const std::vector<Index>& GetOrigVarIdx() const { return orig_var_idx_; }

   private:
    MilpModel& model_;

    Index orig_n_vars_;
    std::vector<Index> orig_var_idx_;

    BitMask con_mask_;
    BitMask var_mask_;
    std::vector<Index> deleted_cons_;
    std::vector<Index> deleted_vars_;

    std::vector<Bounds> activities_;
};

}  // namespace reshala

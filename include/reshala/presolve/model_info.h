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
    ModelInfo(MilpModel& model)
        : model_(model), con_mask(model.GetNCons()), var_mask(model.GetNVars()) {}

    const MilpModel& GetModel() const { return model_; }
    MilpModel& GetModel() { return model_; }

    void MaskCon(Index ic) {
        con_mask.Set(ic);
        deleted_cons.push_back(ic);
    }
    void MaskVar(Index iv) {
        var_mask.Set(iv);
        deleted_vars.push_back(iv);
    }
    inline bool GetConMask(Index ic) const { return con_mask.Get(ic); }
    inline bool GetVarMask(Index ic) const { return var_mask.Get(ic); }

    void CompressVars();

   private:
    MilpModel& model_;

    BitMask con_mask;
    BitMask var_mask;
    std::vector<Index> deleted_cons;
    std::vector<Index> deleted_vars;

    void CalcActivities();
};

}  // namespace reshala

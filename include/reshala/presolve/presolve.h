#pragma once

#include <memory>

#include "reshala/model/milp_model.h"
#include "reshala/presolve/transforms.h"
#include "reshala/presolve/rules.h"

namespace reshala {

enum class ReductionType {
    kFixed,
};

class PostsolveData {
   public:
    void Push(ReductionType type, const std::vector<Index>& indices_,
              const std::vector<Scalar> scalars) {}
    void GetCurrent(int x);

   private:
    std::vector<ReductionType> types;
    std::vector<Index> indices;
    std::vector<Scalar> scalars;
};

class Presolver {
   public:
    Presolver(MilpModel& model);

    void Presolve();
    Solution Postsolve(const Solution&);

   private:
    MilpModel& model_;

    std::vector<std::unique_ptr<Rule>> rules_;
    std::vector<std::unique_ptr<Transform>> transforms_;

    std::vector<Bounds> activities_;

    PostsolveData data_;
};

}  // namespace reshala

#pragma once

#include "reshala/reshala.h"

namespace reshala {

class LpWriter {
   public:
    LpWriter(Reshala& ctx);

    void PrintModel(std::ostream& os, const MilpModel& model) const;

   private:
    Reshala& ctx_;

    void PrintObj(std::ostream& os, const Objective& obj) const;
    void PrintCon(std::ostream& os, const SparseVector& lhs, const std::string& exp,
                  Scalar rhs) const;
};

}  // namespace reshala

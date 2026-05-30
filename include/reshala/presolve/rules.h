#pragma once

#include <memory>
#include <string>

#include "reshala/presolve/model_tracker.h"
#include "reshala/presolve/utils.h"

namespace reshala {

class Rule {
   public:
    Rule(RuleType t) : type(t) {}
    virtual ~Rule() = default;

    virtual RuleResult Apply(ModelTracker& tracker) = 0;

    virtual std::string Name() const = 0;
    RuleType type;
};

class Rule31 : public Rule {
   public:
    Rule31(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "3.1 RedCon"; }
};

class Rule32 : public Rule {
   public:
    Rule32(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "3.2 BndStr"; }
};

class Rule33 : public Rule {
   public:
    Rule33(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "3.3 CoefStr"; }
};

class Rule35 : public Rule {
   public:
    Rule35(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "3.5 Scaling"; }
};

class Rule36 : public Rule {
   public:
    Rule36(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "3.6 SimProb"; }
};

class Rule41 : public Rule {
   public:
    Rule41(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "4.1 FixVar"; }
};

class Rule44 : public Rule {
   public:
    Rule44(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "4.4 DualFix"; }
};

class Rule52 : public Rule {
   public:
    Rule52(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "5.2 ParRows"; }

   private:
    std::vector<Scalar> bin_scales;
    inline Index HashRow(const ModelTracker& tracker, Index ic) const;
    bool Parallel(const SparseVector& sv1, Scalar scale1, const SparseVector& sv2,
                  Scalar scale2) const;
};

class Rule72 : public Rule {
   public:
    Rule72(RuleType t) : Rule(t) {}
    RuleResult Apply(ModelTracker& tracker);
    std::string Name() const { return "7.2 Probing"; }
};

}  // namespace reshala

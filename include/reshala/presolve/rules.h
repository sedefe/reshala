#pragma once

#include <memory>
#include <string>

#include "reshala/presolve/model_tracker.h"
#include "reshala/presolve/utils.h"

namespace reshala {

class Rule {
   public:
    Rule(const std::string& name, RuleType t) : name_(name), type(t) {}
    virtual ~Rule() = default;
    const std::string& GetName() const { return name_; }

    virtual RuleResult Apply(ModelTracker& tracker) = 0;

    RuleType type;

   protected:
    const std::string name_;
};

class Rule31 : public Rule {
   public:
    Rule31(RuleType t) : Rule("3.1 RedCon", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule32 : public Rule {
   public:
    Rule32(RuleType t) : Rule("3.2 BndStr", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule33 : public Rule {
   public:
    Rule33(RuleType t) : Rule("3.3 CoefStr", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule35 : public Rule {
   public:
    Rule35(RuleType t) : Rule("3.5 Scaling", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule36 : public Rule {
   public:
    Rule36(RuleType t) : Rule("3.6 SimProb", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule41 : public Rule {
   public:
    Rule41(RuleType t) : Rule("4.1 FixVar", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule44 : public Rule {
   public:
    Rule44(RuleType t) : Rule("4.4 DualFix", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule46 : public Rule {
   public:
    Rule46(RuleType t) : Rule("4.6 SimSub", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule47 : public Rule {
   public:
    Rule47(RuleType t) : Rule("4.7 RmSlack", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule48 : public Rule {
   public:
    Rule48(RuleType t) : Rule("4.8 Int2Bin", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule52 : public Rule {
   public:
    Rule52(RuleType t) : Rule("5.2 ParRows", t) {}
    RuleResult Apply(ModelTracker& tracker);

   private:
    std::vector<Scalar> bin_scales;
    inline Index HashRow(const ModelTracker& tracker, Index ic) const;
    bool Parallel(const SparseVector& sv1, Scalar scale1, const SparseVector& sv2,
                  Scalar scale2) const;
};

class Rule72 : public Rule {
   public:
    Rule72(RuleType t) : Rule("7.2 Probing", t) {}
    RuleResult Apply(ModelTracker& tracker);
};

}  // namespace reshala

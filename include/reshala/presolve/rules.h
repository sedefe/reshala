#pragma once

#include <memory>
#include <string>

#include "reshala/presolve/model_tracker.h"
#include "reshala/presolve/utils.h"

namespace reshala {

class Rule {
   public:
    Rule(const std::string& name) : name_(name) {}
    virtual ~Rule() = default;
    const std::string& GetName() const { return name_; }

   protected:
    const std::string name_;
};

class RowRule : public Rule {
   public:
    RowRule(const std::string& name) : Rule(name) {}
    virtual RuleResult Apply(ModelTracker& tracker, Index ic) = 0;
};

class ColRule : public Rule {
   public:
    ColRule(const std::string& name) : Rule(name) {}
    virtual RuleResult Apply(ModelTracker& tracker, Index iv) = 0;
};

class HardRule : public Rule {
   public:
    HardRule(const std::string& name) : Rule(name) {}
    virtual RuleResult Apply(ModelTracker& tracker) = 0;
};

class Rule31 : public RowRule {
   public:
    Rule31() : RowRule("3.1 RedCon") {}
    RuleResult Apply(ModelTracker& tracker, Index ic);
};

class Rule32 : public RowRule {
   public:
    Rule32() : RowRule("3.2 BndStr") {}
    RuleResult Apply(ModelTracker& tracker, Index ic);
};

class Rule33 : public RowRule {
   public:
    Rule33() : RowRule("3.3 CoefStr") {}
    RuleResult Apply(ModelTracker& tracker, Index ic);
};

class Rule35 : public RowRule {
   public:
    Rule35() : RowRule("3.5 Scaling") {}
    RuleResult Apply(ModelTracker& tracker, Index ic);
};

class Rule36 : public RowRule {
   public:
    Rule36() : RowRule("3.6 SimProb") {}
    RuleResult Apply(ModelTracker& tracker, Index ic);
};

class Rule41 : public ColRule {
   public:
    Rule41() : ColRule("4.1 FixVar") {}
    RuleResult Apply(ModelTracker& tracker, Index iv);
};

class Rule44 : public ColRule {
   public:
    Rule44() : ColRule("4.4 DualFix") {}
    RuleResult Apply(ModelTracker& tracker, Index iv);
};

class Rule46 : public RowRule {
   public:
    Rule46() : RowRule("4.6 SimSub") {}
    RuleResult Apply(ModelTracker& tracker, Index ic);
};

class Rule47 : public HardRule {
   public:
    Rule47() : HardRule("4.7 RmSlack") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule48 : public ColRule {
   public:
    Rule48() : ColRule("4.8 Int2Bin") {}
    RuleResult Apply(ModelTracker& tracker, Index iv);
};

class Rule52 : public HardRule {
   public:
    Rule52() : HardRule("5.2 ParRows") {}
    RuleResult Apply(ModelTracker& tracker);

   private:
    std::vector<Scalar> bin_scales;
    inline Index HashRow(const ModelTracker& tracker, Index ic) const;
    bool Parallel(const SparseVector& sv1, Scalar scale1, const SparseVector& sv2,
                  Scalar scale2) const;
};

class Rule53 : public HardRule {
    const size_t kMaxNzs = 70;
    bool first_run = true;

   public:
    Rule53() : HardRule("5.3 NzCancel") {}
    RuleResult Apply(ModelTracker& tracker);

   private:
    struct Key {
        Index iv1;
        Index iv2;
        Scalar ratio;
        bool operator==(const Key& o) const {
            return iv1 == o.iv1 && iv2 == o.iv2 && ratio == o.ratio;
        }
    };
    struct KeyHash {
        std::size_t operator()(const Key& k) const noexcept {
            std::size_t h = std::hash<Index>{}(k.iv1);
            h ^= std::hash<Index>{}(k.iv2) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h ^= std::hash<Scalar>{}(k.ratio) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            return h;
        }
    };
    struct Value {
        Index ic;
        Scalar a1;
        Scalar a2;
    };
    std::unordered_map<Key, Value, KeyHash> eq_hash_map;

    void FillEqHashMap(const ModelTracker& tracker);
    Index PairSearch(ModelTracker& tracker);
};

class Rule72 : public HardRule {
   public:
    Rule72() : HardRule("7.2 Probing") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule76 : public HardRule {
   public:
    Rule76() : HardRule("7.6 ImplInt") {}
    RuleResult Apply(ModelTracker& tracker);
};

}  // namespace reshala

#pragma once

#include <memory>
#include <string>
#include <unordered_set>

#include "reshala/presolve/model_tracker.h"
#include "reshala/presolve/utils.h"

namespace reshala {

class Rule {
   public:
    Rule(const std::string& name) : name_(name) {}
    virtual ~Rule() = default;
    const std::string& GetName() const { return name_; }

    virtual RuleResult Apply(ModelTracker& tracker) = 0;

   protected:
    const std::string name_;
};

class Rule31 : public Rule {
   public:
    Rule31() : Rule("3.1 RedCon") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule32 : public Rule {
   public:
    Rule32() : Rule("3.2 BndStr") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule33 : public Rule {
   public:
    Rule33() : Rule("3.3 CoefStr") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule35 : public Rule {
   public:
    Rule35() : Rule("3.5 Scaling") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule36 : public Rule {
   public:
    Rule36() : Rule("3.6 SimProb") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule41 : public Rule {
   public:
    Rule41() : Rule("4.1 FixVar") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule44 : public Rule {
   public:
    Rule44() : Rule("4.4 DualFix") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule46 : public Rule {
   public:
    Rule46() : Rule("4.6 SimSub") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule47 : public Rule {
   public:
    Rule47() : Rule("4.7 RmSlack") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule48 : public Rule {
   public:
    Rule48() : Rule("4.8 Int2Bin") {}
    RuleResult Apply(ModelTracker& tracker);
};

class Rule52 : public Rule {
   public:
    Rule52() : Rule("5.2 ParRows") {}
    RuleResult Apply(ModelTracker& tracker);

   private:
    std::vector<Scalar> bin_scales;
    inline Index HashRow(const ModelTracker& tracker, Index ic) const;
    bool Parallel(const SparseVector& sv1, Scalar scale1, const SparseVector& sv2,
                  Scalar scale2) const;
};

class Rule53 : public Rule {
    const size_t kMaxNzs = 70;
    bool first_run = true;

   public:
    Rule53() : Rule("5.3 NzCancel") {}
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

class Rule72 : public Rule {
    const Index kBudget = 1000;

   public:
    Rule72() : Rule("7.2 Probing") {}
    RuleResult Apply(ModelTracker& tracker);

   private:
    struct BinValue {
        Index iv;
        bool val;
        bool operator==(const BinValue& o) const { return iv == o.iv && val == o.val; }
        bool operator<(const BinValue& o) const {
            if (iv != o.iv) return iv < o.iv;
            return val < o.val;
        }
    };
    struct BinValueHash {
        std::size_t operator()(const BinValue& bv) const noexcept {
            std::size_t h = std::hash<Index>{}(bv.iv);
            if (bv.val) h ^= 0x3AB0D;
            return h;
        }
    };

    std::unordered_map<BinValue, std::unordered_set<BinValue, BinValueHash>, BinValueHash>
        cliq_implics;
    std::vector<std::array<Index, 4>> scores;

    Index FindCliques(ModelTracker& tracker);
    void Score(ModelTracker& tracker);
};

class Rule76 : public Rule {
   public:
    Rule76() : Rule("7.6 ImplInt") {}
    RuleResult Apply(ModelTracker& tracker);
};

}  // namespace reshala

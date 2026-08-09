#include "reshala/heuristics/utils.h"

namespace reshala {

std::string FixingType2Str(FixingType type) {
    switch (type) {
        case FixingType::kAll:
            return "All";
        case FixingType::kInts:
            return "Ints";
        case FixingType::kNone:
            return "None";
        default:
            assert(false && "Unknows rounding type");
            return "";
    }
}

void Fixing(FixingType type, MilpModel &model, const std::vector<Scalar> &relaxed_x) {
    Index n = model.GetNVars();

    Index n_fixed_integers = 0;
    Index n_fixed_continuous = 0;

    if (type != FixingType::kNone) {
        for (Index iv = 0; iv < n; ++iv) {
            if (model.GetIntegrality(iv)) {
                if (MinFraction(relaxed_x[iv]) <= kEpsZero) {
                    auto round_x = Round(relaxed_x[iv]);
                    model.SetBounds(iv, {round_x, round_x});
                    n_fixed_integers++;
                }
            } else {
                if (type == FixingType::kAll) {
                    auto bnd = model.GetBounds(iv);
                    if (WeakLe(relaxed_x[iv], bnd.le)) {
                        model.SetBounds(iv, {bnd.le, bnd.le});
                        n_fixed_continuous++;
                    } else if (WeakGe(relaxed_x[iv], bnd.ri)) {
                        model.SetBounds(iv, {bnd.ri, bnd.ri});
                        n_fixed_continuous++;
                    }
                }
            }
        }
    }

    std::cout << "Fixed " << n_fixed_integers << " integers and " << n_fixed_continuous
              << " continuous\n";
}

}  // namespace reshala

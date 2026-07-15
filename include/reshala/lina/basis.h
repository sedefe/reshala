#pragma once

#include <vector>

#include "reshala/types.h"

namespace reshala {

class LpBasis {
   public:
    LpBasis() {}
    LpBasis(Index m, Index n) {
        basis.resize(m);
        non_basis.resize(n);
        index2nb.resize(m + n, -1);
        Reset();
    }
    void Reset() {
        Index m = basis.size();
        Index n = non_basis.size();
        for (Index ic = 0; ic < m; ic++) {
            basis[ic] = n + ic;
        }
        index2nb.assign(n + m, -1);
        for (Index iv = 0; iv < n; iv++) {
            non_basis[iv] = iv;
            index2nb[iv] = iv;
        }
    }

    inline void Swap(Index iv_leaving, Index iv_entering) {
        index2nb[basis[iv_leaving]] = iv_entering;
        index2nb[non_basis[iv_entering]] = -1;
        std::swap(basis[iv_leaving], non_basis[iv_entering]);
    }

    inline const std::vector<Index>& Basis() const { return basis; }
    inline const std::vector<Index>& NonBasis() const { return non_basis; }
    inline const std::vector<Index>& Index2Nb() const { return index2nb; }

   private:
    std::vector<Index> basis;
    std::vector<Index> non_basis;
    std::vector<Index> index2nb;
};

std::ostream& operator<<(std::ostream& os, const LpBasis& basis);

}  // namespace reshala

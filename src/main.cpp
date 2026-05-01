#include <cstdio>
#include <iostream>

#include "reshala/linalg/sparse_matrix.h"

using namespace reshala;

int main(int argc, char** argv) {
    Index m = 3;
    Index n = 3;
    SparseColMatrix scm(m, n);

    std::vector<SparseVector>& cols = scm.getCols();
    cols[0].push(0, 1);
    cols[1].push(1, 1);
    cols[2].push(2, 1);

    SparseVector sv(n, {0, 1, 2}, {2, 3, 4});
    DenseVector res(n);

    MulScmSv(scm, sv, res);

    std::cout << sv;
    std::cout << res;
}

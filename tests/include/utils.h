#pragma once

#include <cmath>
#include <fstream>

#include "reshala/types.h"

namespace reshala {

class CoutSuppressor {
   private:
    std::streambuf* original_cout;
    std::ofstream null_stream;

   public:
    CoutSuppressor() {
        original_cout = std::cout.rdbuf();
        null_stream.open("/dev/null");
        std::cout.rdbuf(null_stream.rdbuf());
    }

    ~CoutSuppressor() { std::cout.rdbuf(original_cout); }
};

bool CompareScalars(Scalar a, Scalar b) {
    if (std::isnan(a) && std::isnan(b)) {
        return true;
    }
    if (std::isinf(a) && std::isinf(b) && ((a > 0 && b > 0) || (a < 0 && b < 0))) {
        return true;
    }

    Scalar rel_tol = 1e-4;
    auto diff = std::fabs(a - b);
    auto max_val = std::fmax(std::fabs(a), std::fabs(b));

    if (max_val <= rel_tol) {
        return diff <= rel_tol;
    }

    return diff / max_val <= rel_tol;
}

}  // namespace reshala

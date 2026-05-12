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

bool CompareScalars(Scalar x_desired, Scalar x_real) {
    if (std::isnan(x_desired) && std::isnan(x_real)) {
        return true;
    }
    if (std::isinf(x_desired) && std::isinf(x_real) && ((x_desired > 0 && x_real > 0) || (x_desired < 0 && x_real < 0))) {
        return true;
    }

    Scalar rel_tol = 1e-4;
    auto diff = std::fabs(x_real- x_desired);
    auto den = std::fmax(std::fabs(x_desired), rel_tol);

    return diff / den <= rel_tol;
}

}  // namespace reshala

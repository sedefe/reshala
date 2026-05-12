#pragma once

#include <iomanip>

#define FMT(w, p) std::setw(w) << std::setprecision(p) << std::defaultfloat
#define FMT_DEFAULT std::defaultfloat

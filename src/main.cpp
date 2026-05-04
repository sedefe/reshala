#include <cstdio>
#include <iostream>

#include "reshala/io/lp/lp_reader.h"
#include "reshala/lp/dual_simplex.h"

using namespace reshala;

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Usage: %s FILE\n", argv[0]);
        exit(0);
    }

    LpReader reader;
    reader.read(argv[1]);
    auto model = reader.GetModel();

    std::cout << model;

    DualRevisedSimplex drs(model);
}

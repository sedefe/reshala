#include "reshala/io/reader.h"
#include "reshala/lp/dual_simplex.h"
#include "reshala/presolve/presolve.h"

using namespace reshala;

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Usage: %s FILE\n", argv[0]);
        exit(0);
    }

    MilpModel model;
    auto read_status = Read(argv[1], model);
    printf("Reading %s: %s\n", argv[1], FileReadResult2Str(read_status).c_str());
    // std::cout << model;

    Presolver presolver(model);
    presolver.Presolve();

    DualSimplex drs(model);
    const Solution& sol = drs.Solve();
    printf("Status: %s, Obj: %.2f\n", LpStatus2Str(sol.status).c_str(), sol.y);
}

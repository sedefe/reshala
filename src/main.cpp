#include "reshala/io/io.h"
#include "reshala/lp/dual_simplex.h"
#include "reshala/presolve/presolve.h"

using namespace reshala;

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Usage: %s FILE\n", argv[0]);
        exit(0);
    }

    Io io;
    auto read_status = io.Read(argv[1]);
    MilpModel& model = io.GetModel();
    printf("Reading %s: %s\n", argv[1], FileReadResult2Str(read_status).c_str());
    if (read_status != FileReadStatus::kOk) {
        exit(0);
    }
    // std::cout << model;

    Presolver presolver(model);
    presolver.Presolve();

    DualSimplex drs(model);
    const Solution& sol = drs.Solve();
    printf("Status: %s, Obj: %.2f\n", LpStatus2Str(sol.status).c_str(), sol.y);
    // for (Index i = 0; i < model.GetNVars(); i++) {
    //     printf("x[%2d] = %5.2f\n", i, sol.x[i]);
    // }
}

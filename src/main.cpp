#include "reshala/io/lp/lp_reader.h"
#include "reshala/io/mps/mps_reader.h"
#include "reshala/lp/dual_simplex.h"
#include "reshala/presolve/presolve.h"

using namespace reshala;

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Usage: %s FILE\n", argv[0]);
        exit(0);
    }

    // LpReader reader;
    MpsReader reader;
    auto read_status = reader.Read(argv[1]);
    printf("Reading %s: %s\n", argv[1], FileReadResult2Str(read_status).c_str());
    auto model = reader.GetModel();
    std::cout << model;
    exit(0);

    Presolver presolver(model);
    presolver.Presolve();

    DualSimplex drs(model);
    const Solution& sol = drs.Solve();
    printf("Status: %s, Obj: %.2f\n", LpStatus2Str(sol.status).c_str(), sol.y);
}

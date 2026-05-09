#include "reshala/io/io.h"
#include "reshala/milp/milp.h"

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
    const MilpModel model_copy = model;
    // std::cout << model;
    printf("%d vars, %d cons\n", model.GetNVars(), model.GetNCons());

    MilpSolver milp(model);

    auto start = std::chrono::high_resolution_clock::now();
    const Solution& sol = milp.Solve();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    printf("Solved in %.3f ms\n", duration.count() / 1e3);

    printf("Status: %s, Obj: %8.5g\n", LpStatus2Str(sol.status).c_str(), sol.y);
    if (sol.status == LpStatus::kOptimal) {
        auto rep = model_copy.GetFeasReport(sol.x);
        printf("Violations: int %.5f, bnd %.5f, con %.5f\n", rep.max_int_infeas, rep.max_bnd_infeas,
               rep.max_con_infeas);
        io.PrintValues(std::cout, sol.x);
    }
}

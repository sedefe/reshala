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
    // std::cout << model;

    printf("%d vars, %d cons\n", model.GetNVars(), model.GetNCons());

    MilpSolver milp(model);

    auto start = std::chrono::high_resolution_clock::now();
    const Solution& sol = milp.Solve();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    printf("Solved in %.3f ms\n", duration.count() / 1e3);

    printf("Status: %s, Obj: %.2f\n", LpStatus2Str(sol.status).c_str(), sol.y);
    if (sol.status == LpStatus::kOptimal) {
        io.PrintValues(std::cout, sol.x);
    }
}

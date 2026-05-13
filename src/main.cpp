#include "reshala/io/io.h"
#include "reshala/milp/milp.h"

using namespace reshala;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " FILE\n";
        exit(0);
    }

    Io io;
    std::cout << "Reading " << argv[1] << "\n";
    auto read_status = io.Read(argv[1]);
    MilpModel& model = io.GetModel();
    if (read_status != FileReadStatus::kOk) {
        exit(0);
    }
    const MilpModel model_copy = model;
    // std::cout << model;
    std::cout << model.GetNCons() << " cons, " << model.GetNVars() << " vars\n";

    MilpSolver milp(model);

    auto start = std::chrono::high_resolution_clock::now();
    const Solution& sol = milp.Solve();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Solved in " << duration.count() / 1e3 << " ms\n";

    std::cout << "Status: " << LpStatus2Str(sol.status) << ", objective: " << FMT(-10, 5) << sol.y
              << "\n";
    if (sol.status == LpStatus::kOptimal) {
        auto rep = model_copy.GetFeasReport(sol.x);
        std::cout << "Violations: int " << rep.max_int_infeas << ", bnd " << rep.max_bnd_infeas
                  << ", con " << rep.max_con_infeas << "\n";
        io.PrintValues(std::cout, sol.x);
    }
}

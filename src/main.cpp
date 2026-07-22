#include "reshala/io/io.h"
#include "reshala/milp/milp.h"

using namespace reshala;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " FILE\n";
        exit(0);
    }
    std::cout << "Build: " << kBuildType.c_str() << "\n";

    Io io;
    std::cout << "Reading " << argv[1] << "\n";
    auto read_status = io.Read(argv[1]);
    if (read_status != FileReadStatus::kOk) {
        std::cerr << FileReadStatus2Str(read_status) << "\n";
        exit(0);
    }
    MilpModel& model = io.GetModel();
    const MilpModel model_copy = model;
    // std::cout << model;
    std::cout << model.StatString() << "\n";

    MilpSolver milp(model);
    auto [sol, t_solve] = MEASURE_TIME(milp.Solve());
    std::cout << "Solved in " << t_solve << " ms\n";

    std::cout << "Status: " << LpStatus2Str(sol.status) << ", objective: " << FMT(-10, 5) << sol.y
              << "\n";
    milp.PrintStats(std::cout);
    if (sol.status == LpStatus::kOptimal) {
        std::cout << "=== Checking ===\n";
        auto y = model_copy.GetObj().evaluate(sol.x);
        std::cout << "Objective: " << FMT(-10, 5) << y << "\n";
        auto rep = model_copy.GetFeasReport(sol.x);
        std::cout << rep << "\n";
        io.PrintValues(std::cout, sol.x);
    }
}

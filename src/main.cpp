#include "reshala/reshala.h"

using namespace reshala;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " FILE\n";
        exit(0);
    }
    std::cout << "Build: " << kBuildType.c_str() << "\n";

    Reshala reshala;
    std::cout << "Reading " << argv[1] << "\n";

    auto read_status = reshala.Read(argv[1]);
    if (read_status != FileReadStatus::kOk) {
        std::cerr << FileReadStatus2Str(read_status) << "\n";
        exit(0);
    }

    const MilpModel model = reshala.GetModel();  // For verification
    std::cout << model.StatString() << "\n";

    auto [sol, t_solve] = MEASURE_TIME(reshala.Solve());
    std::cout << "Solved in " << t_solve << " ms\n";

    std::cout << "Status: " << LpStatus2Str(sol.status) << ", objective: " << FMT(-10, 5) << sol.y
              << "\n";
    reshala.PrintStats(std::cout);
    if (sol.status == LpStatus::kOptimal) {
        std::cout << "=== Checking ===\n";
        auto y = model.GetObj().evaluate(sol.x);
        std::cout << "Objective: " << FMT(-10, 5) << y << "\n";
        auto rep = model.GetFeasReport(sol.x);
        std::cout << rep << "\n";
        // io.PrintValues(std::cout, sol.x);
    }
}

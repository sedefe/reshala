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
    if (read_status != FileReadStatus::kOk) {
        std::cerr << FileReadStatus2Str(read_status) << "\n";
        exit(0);
    }
    MilpModel& model = io.GetModel();
    const MilpModel model_copy = model;
    // std::cout << model;
    std::cout << model.GetNCons() << " cons, " << model.GetNVars() << " vars\n";

    MilpSolver milp(model);
    auto [sol, duration] = MEASURE_TIME(milp.Solve());
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

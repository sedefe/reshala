#include <chrono>

#include "reshala/io/io.h"
#include "reshala/milp/milp.h"
#include "utils.h"

using namespace reshala;

struct TestCase {
    std::string name;
    std::string status;
    Scalar y;

    Solution sol;
    std::chrono::microseconds time;
    FeasibilityReport report;
};

std::vector<TestCase> ReadTestCases(const std::string& csv_path) {
    std::vector<TestCase> tests;
    std::ifstream file(csv_path);

    if (!file.is_open()) {
        std::cerr << "Failed to open: " << csv_path << std::endl;
        return tests;
    }

    std::string line;

    std::getline(file, line);
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        TestCase tc;
        std::string token;

        std::getline(ss, tc.name, ',');
        std::getline(ss, token, ',');
        tc.status = token;
        std::getline(ss, token, ',');
        if (tc.status == "Optimal") {
            tc.y = std::stod(token);
        } else {
            tc.y = kInf;
        }
        std::getline(ss, token, ',');
        bool is_used = std::stoi(token);
        if (is_used) {
            tests.push_back(tc);
        }
    }

    return tests;
}

void RunTest(TestCase& tc) {
    Io io;
    std::filesystem::path file_path("tests/models/" + tc.name + ".mps");
    [[maybe_unused]] auto read_status = io.Read(file_path.c_str());
    assert(read_status == FileReadStatus::kOk);

    MilpModel& model = io.GetModel();
    MilpModel model_copy = model;
    MilpSolver solver(model);
    Solution sol;

    auto [_, time] = MEASURE_TIME({
        CoutSuppressor suppressor;
        tc.sol = solver.Solve();
    });
    tc.time = time;

    if (tc.sol.status == LpStatus::kOptimal) {
        tc.report = model_copy.GetFeasReport(tc.sol.x);
    }
}

int main() {
#if defined(NDEBUG)
    printf("Build: release\n");
#elif defined(DEBUG)
    printf("Build: debug\n");
#else
    printf("Build: X3\n");
#endif

    std::vector<TestCase> test_cases = ReadTestCases("tests/models/results.csv");

    size_t n_tests = test_cases.size();
    size_t n_good = 0;

    for (auto& tc : test_cases) {
        bool is_good = true;
        printf("%-20s: ", tc.name.c_str());
        fflush(stdout);

        RunTest(tc);
        printf("%6.3f sec ", tc.time.count() / 1e6);

        auto status = LpStatus2Str(tc.sol.status);
        printf("%12s ", status.c_str());
        if (status == tc.status) {
            printf("[√] ");
        } else {
            printf("[X] (%12s)", tc.status.c_str());
            is_good = false;
        }

        bool compare_y = CompareScalars(tc.y, tc.sol.y);
        printf("Obj: %10.5g ", tc.sol.y);
        if (compare_y) {
            printf("[√] ");
        } else {
            printf("[X] (%.5g)", tc.y);
            is_good = false;
        }

        if (tc.sol.status == LpStatus::kOptimal) {
            printf("Feasibility (I/B/C): ");
            auto& [aii, abi, aci, rbi, rci] = tc.report;
            for (Scalar* s : {&aii, &abi, &aci}) {
                if (*s <= kEpsZero) {
                    printf("√");
                } else {
                    printf("X");
                    is_good = false;
                }
            }
        }

        printf("\n");

        n_good += is_good;
    }

    printf("%lu/%lu tests passed\n", n_good, n_tests);
}

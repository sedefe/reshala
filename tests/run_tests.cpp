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

bool RunTest(TestCase& tc) {
    Io io;
    std::filesystem::path file_path("tests/models/" + tc.name + ".mps");
    auto read_status = io.Read(file_path.c_str());

    MilpModel& model = io.GetModel();
    MilpModel model_copy = model;
    MilpSolver solver(model);
    Solution sol;

    auto start = std::chrono::high_resolution_clock::now();
    {
        CoutSuppressor suppressor;
        tc.sol = solver.Solve();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    tc.time = time;

    if (tc.sol.status == LpStatus::kOptimal) {
        tc.report = model_copy.GetFeasReport(tc.sol.x);
    }

    return true;
}

int main() {
    std::vector<TestCase> test_cases = ReadTestCases("tests/models/results.csv");

    size_t n_tests = test_cases.size();
    size_t n_good = 0;

    for (auto& tc : test_cases) {
        bool is_good = true;
        printf("%-20s: ", tc.name.c_str());
        fflush(stdout);

        bool res = RunTest(tc);
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
            auto& [iv, bv, cv] = tc.report;
            for (Scalar* p : {&iv, &bv, &cv}) {
                if (iv <= kEpsZero) {
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

#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "reshala/io/io.h"
#include "reshala/milp/milp.h"

using namespace reshala;

class CoutSuppressor {
   private:
    std::streambuf* original_cout;
    std::ofstream null_stream;

   public:
    CoutSuppressor() {
        original_cout = std::cout.rdbuf();
        null_stream.open("/dev/null");
        std::cout.rdbuf(null_stream.rdbuf());
    }

    ~CoutSuppressor() { std::cout.rdbuf(original_cout); }
};

struct TestCase {
    std::string name;
    std::string status;
    Scalar y;
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

Solution RunTest(const std::string& model_name) {
    Io io;
    std::filesystem::path file_path("tests/models/" + model_name + ".mps");
    auto read_status = io.Read(file_path.c_str());
    printf("Model %s:\n", model_name.c_str());

    MilpModel& model = io.GetModel();
    MilpSolver solver(model);
    Solution sol;

    auto start = std::chrono::high_resolution_clock::now();
    {
        CoutSuppressor suppressor;
        sol = solver.Solve();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    printf("Solved in %.3f ms\n", duration.count() / 1e3);

    return sol;
}

int main() {
    std::vector<TestCase> test_cases = ReadTestCases("tests/models/results.csv");

    for (auto& test : test_cases) {
        auto sol = RunTest(test.name);
        printf("%20s: Status: %12s (%12s), Obj: %.5g (%.5g)\n", test.name.c_str(),
               LpStatus2Str(sol.status).c_str(), test.status.c_str(), sol.y, test.y);
    }
}

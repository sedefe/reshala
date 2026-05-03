#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "reshala/io/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

enum class LpReadResult { kOk, kFsError, kParseError };
enum class ParseState { kObj, kCon, kBnd, kBin, kGen, kDon, kNon };

class LpReader {
   public:
    LpReader() {}

    LpReadResult read(const char* fname);
    MilpModel& GetModel() { return model; }

   private:
    MilpModel model;
    NameMapper var_names;

    std::vector<std::string> tokenize_line(const std::string& line);

    void parse_objective(const std::vector<std::string>&);
    void parse_constraint(const std::vector<std::string>&);
    void parse_bounds(const std::vector<std::string>&);
    void parse_binaries(const std::vector<std::string>&);
    void parse_generals(const std::vector<std::string>&);
};

}  // namespace reshala

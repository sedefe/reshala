#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "reshala/io/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

class AbstractReader {
   public:
    AbstractReader(MilpModel& model, Names& names) : model_(model), names_(names) {}
    virtual FileReadStatus Parse(std::ifstream& file) = 0;

   protected:
    MilpModel& model_;
    Index line_number = 0;

    void ThrowParseError(const std::string& message) const {
        throw std::runtime_error("Line " + std::to_string(line_number) + ": " + message + "\n");
    }

    const std::string& GetObjName() const { return names_.obj; }
    void SetObjName(const std::string& name) { names_.obj = name; }
    Index GetConIndex(const std::string& name, bool can_add) {
        auto res = names_.cons.GetIndex(name, can_add);
        if (res < 0) ThrowParseError("Unexpected constraint: " + name);
        return res;
    }
    Index GetVarIndex(const std::string& name, bool can_add) {
        auto res = names_.vars.GetIndex(name, can_add);
        if (res < 0) ThrowParseError("Unexpected variable: " + name);
        return res;
    }

    const Names& GetNames() const { return names_; }
    Index GetNCons() const { return names_.cons.Size(); }
    Index GetNVars() const { return names_.vars.Size(); }

   private:
    Names& names_;  // Derived readers have no direct access
};

}  // namespace reshala

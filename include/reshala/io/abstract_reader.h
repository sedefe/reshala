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

    Index LookupCon(const std::string& name) const {
        if (auto idx = names_.cons.Lookup(name)) return *idx;
        ThrowParseError("Unexpected constraint: " + name);
        return -1;
    }
    Index GetOrAddCon(const std::string& name) { return names_.cons.GetOrAdd(name); }

    Index LookupVar(const std::string& name) const {
        if (auto idx = names_.vars.Lookup(name)) return *idx;
        ThrowParseError("Unexpected variable: " + name);
        return -1;
    }
    Index GetOrAddVar(const std::string& name) { return names_.vars.GetOrAdd(name); }

    const Names& GetNames() const { return names_; }
    Index GetNCons() const { return names_.cons.Size(); }
    Index GetNVars() const { return names_.vars.Size(); }

   private:
    Names& names_;  // Derived readers have no direct access
};

}  // namespace reshala

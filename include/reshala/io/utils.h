#pragma once

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "reshala/linalg/types.h"
#include "reshala/model/domain.h"

namespace reshala {

enum class FileReadStatus { kOk, kFsError, kParseError };
std::string FileReadResult2Str(FileReadStatus status);

enum class ExpType { kGe, kLe, kEq, kNon };
Bounds ExpType2Bounds(ExpType exp_type, Scalar rhs);

std::string to_lowercase(const std::string& s);

std::vector<std::string> tokenize_line(const std::string& line);

struct NameMapper {
    std::unordered_map<std::string, size_t> name_to_index;
    std::vector<std::string> index_to_name;

    size_t get_index(const std::string& name) {
        auto it = name_to_index.find(name);
        if (it != name_to_index.end()) {
            return it->second;
        }

        size_t new_index = index_to_name.size();
        name_to_index[name] = new_index;
        index_to_name.push_back(name);
        return new_index;
    }

    const std::string& get_name(size_t index) const { return index_to_name.at(index); }

    size_t size() const { return index_to_name.size(); }
};

struct Names {
    NameMapper vars;
    NameMapper cons;
    std::string obj;
};

}  // namespace reshala

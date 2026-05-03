#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "reshala/linalg/types.h"

namespace reshala {

std::string to_lowercase(const std::string& s) {
    std::string res = s;
    for (char& c : res) {
        c = std::tolower((unsigned char)c);
    }
    return res;
}

std::vector<std::string> tokenize_line(const std::string& line) {
    std::vector<std::string> tokens;
    std::string token;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (std::isspace(c)) {
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
        } else {
            token += c;
        }
    }
    if (!token.empty()) tokens.push_back(token);
    return tokens;
}

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

}  // namespace reshala

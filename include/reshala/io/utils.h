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

struct Monom {
    Scalar coeff;
    size_t index;
};

std::vector<Monom> parse_to_monomes(const std::vector<std::string>& tokens, NameMapper& names) {
    std::vector<Monom> result;

    Scalar sign = 1.0;
    Scalar coeff = 1.0;

    for (size_t i = 0; i < tokens.size(); ++i) {
        const std::string& token = tokens[i];
        // printf("[%d] %s\n", i, token.c_str());

        if (token == "+") {
            sign = 1.0;
        } else if (token == "-") {
            sign = -1.0;
        } else if (isdigit(token[0]) || token[0] == '.') {
            Scalar num = std::stod(token);
            coeff = sign * num;
        } else {
            size_t idx = names.get_index(token);
            result.push_back({coeff, idx});

            sign = 1.0;
            coeff = 1.0;
        }
    }
    return result;
}

}  // namespace reshala

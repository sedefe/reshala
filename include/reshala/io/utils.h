#pragma once

#include <string>

std::string to_lowercase(const std::string &s) {
    std::string res = s;
    for (char &c : res) {
        c = std::tolower((unsigned char)c);
    }
    return res;
}

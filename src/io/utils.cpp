#include "reshala/io/utils.h"

#include <assert.h>

namespace reshala {

std::string FileReadStatus2Str(FileReadStatus status) {
    switch (status) {
        case FileReadStatus::kOk:
            return "Ok";
        case FileReadStatus::kFsError:
            return "FS error";
        case FileReadStatus::kParseError:
            return "Parsing error";
        default:
            assert(false && "Unknows status");
            return "";
    }
}

Bounds ExpType2Bounds(ExpType exp_type, Scalar rhs) {
    switch (exp_type) {
        case ExpType::kLe:
            return {-kInf, rhs};
            break;
        case ExpType::kGe:
            return {rhs, kInf};
            break;
        case ExpType::kEq:
            return {rhs, rhs};
            break;
        default:
            assert(false && "Unknown expression type");
            return {kNan, kNan};
    }
}

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

}  // namespace reshala

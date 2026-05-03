#include "reshala/io/lp/lp_reader.h"

#include <iostream>

namespace reshala {

std::vector<std::string> LpReader::tokenize_line(const std::string& line) {
    std::vector<std::string> tokens;
    std::string token;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (std::isspace(c)) {
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
        } else if (c == '+' || c == '-' || c == '<' || c == '>' || c == '=') {
            // Handle multi-char operators like <=, >=, ==
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
            std::string op(1, c);
            if (i + 1 < line.size() && (line[i + 1] == '=' || (c == '<' && line[i + 1] == '='))) {
                op += line[++i];
            }
            tokens.push_back(op);
        } else {
            token += c;
        }
    }
    if (!token.empty()) tokens.push_back(token);
    return tokens;
}

LpReadResult LpReader::read(const char* fname) {
    std::filesystem::path file_path(fname);
    std::ifstream file(file_path);
    if (!file.is_open()) {
        // throw std::runtime_error("Failed to open file: " + std::string(fname));
        return LpReadResult::kFsError;
    }

    std::string line;
    ParseState current_state = ParseState::kNon;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '\\') continue;

        std::string lower_line = to_lowercase(line);

        if (lower_line.find("min") != std::string::npos) {
            model.getObj().sense = Sense::kMin;
            current_state = ParseState::kObj;
            continue;
        } else if (lower_line.find("max") != std::string::npos) {
            model.getObj().sense = Sense::kMax;
            current_state = ParseState::kObj;
            continue;
        } else if (lower_line.find("so that") != std::string::npos ||
                   lower_line.find("subject to") != std::string::npos) {
            current_state = ParseState::kCon;
            continue;
        } else if (lower_line.find("bounds") != std::string::npos) {
            current_state = ParseState::kBnd;
            continue;
        } else if (lower_line.find("binaries") != std::string::npos) {
            current_state = ParseState::kBin;
            continue;
        } else if (lower_line.find("generals") != std::string::npos) {
            current_state = ParseState::kGen;
            continue;
        } else if (lower_line.find("end") != std::string::npos) {
            current_state = ParseState::kDon;
            break;
        }

        auto tokens = tokenize_line(line);
        if (tokens.empty()) continue;

        switch (current_state) {
            case ParseState::kObj:
                parse_objective(tokens);
                break;
            case ParseState::kCon:
                parse_constraint(tokens);
                break;
            case ParseState::kBnd:
                parse_bounds(tokens);
                break;
            case ParseState::kBin:
                parse_binaries(tokens);
                break;
            case ParseState::kGen:
                parse_generals(tokens);
                break;
            default:
                break;
        }
    }
    file.close();
    return LpReadResult::kOk;
}

void LpReader::parse_objective(const std::vector<std::string>& tokens) {
    auto monomes = parse_to_monomes(tokens, var_names);
    Objective& obj = model.getObj();
    model.getObj().coefficients.resize(var_names.size());
    obj.c0 = 0;
    for (const auto& m : monomes) {
        obj.coefficients[m.index] = m.coeff;
    }
}
void LpReader::parse_constraint(const std::vector<std::string>& tokens) {
    printf("\nCons:\n");
    for (const auto& str : tokens) {
        std::cout << str << std::endl;
    }
}
void LpReader::parse_bounds(const std::vector<std::string>& tokens) {
    printf("\nBounds:\n");
    for (const auto& str : tokens) {
        std::cout << str << std::endl;
    }
}
void LpReader::parse_binaries(const std::vector<std::string>& tokens) {
    printf("\nBins:\n");
    for (const auto& str : tokens) {
        std::cout << str << std::endl;
    }
}
void LpReader::parse_generals(const std::vector<std::string>& tokens) {
    printf("\nGenerals:\n");
    for (const auto& str : tokens) {
        std::cout << str << std::endl;
    }
}

}  // namespace reshala

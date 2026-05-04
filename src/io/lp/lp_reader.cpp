#include "reshala/io/lp/lp_reader.h"

#include <iostream>

namespace reshala {

LpReadResult LpReader::Read(const char* fname) {
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
            model.GetObj().orig_sense = Sense::kMin;
            current_state = ParseState::kObj;
            continue;
        } else if (lower_line.find("max") != std::string::npos) {
            model.GetObj().orig_sense = Sense::kMax;
            current_state = ParseState::kObj;
            continue;
        } else if (lower_line.find("s.t.") != std::string::npos ||
                   lower_line.find("subject to") != std::string::npos) {
            current_state = ParseState::kCon;
            continue;
        } else if (lower_line.find("bounds") != std::string::npos) {
            // Now the size is known
            size_t n_cons = model.GetAr().GetRows().size();
            size_t n_vars = var_names.size();
            model.Resize(n_cons, n_vars);
            model.FinalizeAc();

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
                ParseObjective(tokens);
                break;
            case ParseState::kCon:
                ParseConstraint(tokens);
                break;
            case ParseState::kBnd:
                ParseBounds(tokens);
                break;
            case ParseState::kBin:
                ParseBinaries(tokens);
                break;
            case ParseState::kGen:
                ParseGenerals(tokens);
                break;
            default:
                break;
        }
    }
    file.close();
    return LpReadResult::kOk;
}

void LpReader::ParseObjective(const std::vector<std::string>& tokens) {
    std::vector<Monom> lhs;
    ParseLincomb(tokens, lhs);

    Objective& obj = model.GetObj();
    model.GetObj().coefficients.resize(var_names.size());
    obj.c0 = 0;  // Todo: parse c0
    if (model.GetObj().orig_sense == Sense::kMin) {
        for (const auto& m : lhs) {
            obj.coefficients[m.index] = m.coeff;
        }
    } else {
        for (const auto& m : lhs) {
            obj.coefficients[m.index] = -m.coeff;
        }
    }
}

void LpReader::ParseConstraint(const std::vector<std::string>& tokens) {
    std::vector<Monom> lhs;
    Bounds rhs;
    ParseLincomb(tokens, lhs, tokens.size() - 2);

    const std::string& exp_token = tokens[tokens.size() - 2];
    const std::string& rhs_token = tokens[tokens.size() - 1];

    assert(exp_token[0] == '<' or exp_token[0] == '>' or exp_token[0] == '=');
    assert((exp_token.size() == 1) or (exp_token.size() == 2 and exp_token[1] == '='));
    ExpType type = char2exptype(exp_token[0]);
    Scalar coeff = std::stod(rhs_token);

    switch (type) {
        case ExpType::kLe:
            rhs = {-kInf, coeff};
            break;
        case ExpType::kGe:
            rhs = {coeff, kInf};
            break;
        case ExpType::kEq:
            rhs = {coeff, coeff};
            break;
        default:
            assert(false);
    }

    SparseVector sv(var_names.size());
    sv.Reserve(lhs.size());
    for (const auto& m : lhs) {
        sv.Push(m.index, m.coeff);
    }
    model.PrepareConstraint(sv, rhs);
}

void LpReader::ParseBounds(const std::vector<std::string>& tokens) {
    size_t n = tokens.size();
    assert(n % 3 == 0);
    for (Index i = 0; i < n; i += 3) {
        Index index = var_names.get_index(tokens[i]);
        ExpType type = char2exptype(tokens[i + 1][0]);
        Scalar rhs = std::stod(tokens[i + 2]);
        Bounds& bnd = model.GetBounds(index);
        switch (type) {
            case ExpType::kGe:
                bnd = BoundsIntersection(bnd, {rhs, kInf});
                break;
            case ExpType::kLe:
                bnd = BoundsIntersection(bnd, {-kInf, rhs});
                break;
            case ExpType::kEq:
                bnd = BoundsIntersection(bnd, {rhs, rhs});
                break;
        }
    }
}

void LpReader::ParseBinaries(const std::vector<std::string>& tokens) {
    for (const auto& str : tokens) {
        auto index = var_names.get_index(str);
        model.GetVars().integrality[index] = true;
        Bounds& bnd = model.GetBounds(index);
        bnd = BoundsIntersection(bnd, {0, 1});
    }
}

void LpReader::ParseGenerals(const std::vector<std::string>& tokens) {
    for (const auto& str : tokens) {
        auto index = var_names.get_index(str);
        model.GetVars().integrality[index] = true;
    }
}

void LpReader::ParseLincomb(const std::vector<std::string>& tokens, std::vector<Monom>& lhs,
                            size_t n_tokens) {
    if (n_tokens == -1) n_tokens = tokens.size();
    Scalar sign = 1.0;
    Scalar coeff = 1.0;

    Index i;
    for (i = 0; i < n_tokens; ++i) {
        const std::string& token = tokens[i];

        if (token == "+") {
            assert(token.size() == 1);
            sign = 1.0;
        } else if (token == "-") {
            assert(token.size() == 1);
            sign = -1.0;

        } else if (isdigit(token[0]) || token[0] == '.') {
            coeff = std::stod(token);
        } else {
            size_t idx = var_names.get_index(token);
            lhs.push_back({sign * coeff, idx});

            sign = 1.0;
            coeff = 1.0;
        }
    }
}

}  // namespace reshala

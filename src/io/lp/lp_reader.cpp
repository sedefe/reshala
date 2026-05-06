#include "reshala/io/lp/lp_reader.h"

#include <iostream>

namespace reshala {

FileReadStatus LpReader::Read() {
    std::ifstream file(path_);
    if (!file.is_open()) {
        return FileReadStatus::kFsError;
    }

    std::string line;
    auto current_state = LpParseState::kNon;

    while (std::getline(file, line)) {
        line_number++;
        while (line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '\\') continue;

        auto tokens = tokenize_line(line);
        if (tokens.empty()) continue;
        std::string lower_line = to_lowercase(line);

        if (lower_line == "min" or lower_line == "minimize") {
            model_.GetObj().orig_sense = Sense::kMin;
            current_state = LpParseState::kObj;
            continue;
        } else if (lower_line == "max" or lower_line == "maximize") {
            model_.GetObj().orig_sense = Sense::kMax;
            current_state = LpParseState::kObj;
            continue;
        } else if (lower_line == "s.t." or lower_line == "subject to") {
            current_state = LpParseState::kCon;
            continue;
        } else if (lower_line == "bounds") {
            // Now the size is known
            size_t n_cons = model_.GetAr().GetRows().size();
            size_t n_vars = var_names.size();
            model_.Resize(n_cons, n_vars);
            model_.FinalizeAc();

            current_state = LpParseState::kBnd;
            continue;
        } else if (lower_line == "binaries") {
            current_state = LpParseState::kBin;
            continue;
        } else if (lower_line == "generals") {
            current_state = LpParseState::kGen;
            continue;
        } else if (lower_line == "end") {
            current_state = LpParseState::kDon;
            break;
        }

        switch (current_state) {
            case LpParseState::kObj:
                ParseObjective(tokens);
                break;
            case LpParseState::kCon:
                ParseConstraint(tokens);
                break;
            case LpParseState::kBnd:
                ParseBounds(tokens);
                break;
            case LpParseState::kBin:
                ParseBinaries(tokens);
                break;
            case LpParseState::kGen:
                ParseGenerals(tokens);
                break;
            default:
                break;
        }
    }
    file.close();
    return FileReadStatus::kOk;
}

void LpReader::ParseObjective(const std::vector<std::string>& tokens) {
    std::vector<Monom> lhs;
    ParseLincomb(tokens, lhs);

    Objective& obj = model_.GetObj();
    model_.GetObj().coefficients.resize(var_names.size());
    obj.c0 = 0;  // Todo: parse c0
    if (model_.GetObj().orig_sense == Sense::kMin) {
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

    ParseLincomb(tokens, lhs, tokens.size() - 2);

    const std::string& exp_token = tokens[tokens.size() - 2];
    const std::string& rhs_token = tokens[tokens.size() - 1];

    if (!(exp_token[0] == '<' or exp_token[0] == '>' or exp_token[0] == '=') or
        !((exp_token.size() == 1) or (exp_token.size() == 2 and exp_token[1] == '='))) {
        ThrowParseError("Unsupported expression type: " + exp_token);
    }
    ExpType exp_type = LpChar2ExpType(exp_token[0]);
    Scalar coeff = std::stod(rhs_token);
    Bounds rhs = ExpType2Bounds(exp_type, coeff);

    SparseVector sv(var_names.size());
    sv.Reserve(lhs.size());
    for (const auto& m : lhs) {
        sv.Push(m.index, m.coeff);
    }
    model_.PrepareConstraint(sv, rhs);
}

void LpReader::ParseBounds(const std::vector<std::string>& tokens) {
    size_t n = tokens.size();
    if (n % 3 != 0) {
        ThrowParseError("Can't parse bounds");
    }
    for (Index i = 0; i < n; i += 3) {
        Index index = var_names.get_index(tokens[i]);
        ExpType type = LpChar2ExpType(tokens[i + 1][0]);
        Scalar rhs = std::stod(tokens[i + 2]);
        Bounds& bnd = model_.GetBounds(index);
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
        model_.GetVars().integrality[index] = true;
        Bounds& bnd = model_.GetBounds(index);
        bnd = BoundsIntersection(bnd, {0, 1});
    }
}

void LpReader::ParseGenerals(const std::vector<std::string>& tokens) {
    for (const auto& str : tokens) {
        auto index = var_names.get_index(str);
        model_.GetVars().integrality[index] = true;
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
            sign = 1.0;
        } else if (token == "-") {
            sign = -1.0;
        } else if (isdigit(token[0]) || token[0] == '-' || token[0] == '.') {
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

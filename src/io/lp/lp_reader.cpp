#include "reshala/io/lp/lp_reader.h"

#include <iostream>

namespace reshala {

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
            model.getObj().orig_sense = Sense::kMin;
            current_state = ParseState::kObj;
            continue;
        } else if (lower_line.find("max") != std::string::npos) {
            model.getObj().orig_sense = Sense::kMax;
            current_state = ParseState::kObj;
            continue;
        } else if (lower_line.find("so that") != std::string::npos ||
                   lower_line.find("subject to") != std::string::npos) {
            current_state = ParseState::kCon;
            continue;
        } else if (lower_line.find("bounds") != std::string::npos) {
            finalize(); // Now the size is known
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
    std::vector<Monom> lhs;
    parse_lincomb(tokens, lhs);

    Objective& obj = model.getObj();
    model.getObj().coefficients.resize(var_names.size());
    obj.c0 = 0;  // Todo: parse c0
    if (model.getObj().orig_sense == Sense::kMin) {
        for (const auto& m : lhs) {
            obj.coefficients[m.index] = m.coeff;
        }
    } else {
        for (const auto& m : lhs) {
            obj.coefficients[m.index] = -m.coeff;
        }
    }
}

void LpReader::parse_constraint(const std::vector<std::string>& tokens) {
    std::vector<Monom> lhs;
    Bound rhs;
    parse_lincomb(tokens, lhs, tokens.size() - 2);

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
    sv.reserve(lhs.size());
    for (const auto& m : lhs) {
        sv.push(m.index, m.coeff);
    }
    model.PrepareConstraint(sv, rhs);
}

void LpReader::parse_bounds(const std::vector<std::string>& tokens) {
    size_t n = tokens.size();
    assert(n % 3 == 0);
    for (Index i = 0; i < n; i += 3) {
        Index index = var_names.get_index(tokens[i]);
        ExpType type = char2exptype(tokens[i + 1][0]);
        Scalar rhs = std::stod(tokens[i + 2]);
        Bound& bnd = model.getBounds()[index];
        switch (type) {
            case ExpType::kGe:
                bnd.le = std::max(bnd.le, rhs);
                break;
            case ExpType::kLe:
                bnd.ri = std::min(bnd.ri, rhs);
                break;
            case ExpType::kEq:
                bnd.le = std::max(bnd.le, rhs);
                bnd.ri = std::min(bnd.ri, rhs);
                break;
        }
    }
}

void LpReader::parse_binaries(const std::vector<std::string>& tokens) {
    for (const auto& str : tokens) {
        auto index = var_names.get_index(str);
        model.getIntegrality()[index] = true;
        Bound& bnd = model.getBounds()[index];
        bnd.le = std::max(bnd.le, 0.0);
        bnd.ri = std::min(bnd.ri, 1.0);
    }
}

void LpReader::parse_generals(const std::vector<std::string>& tokens) {
    for (const auto& str : tokens) {
        auto index = var_names.get_index(str);
        model.getIntegrality()[index] = true;
    }
}

void LpReader::finalize() {
    size_t n_cons = model.getAr().getRows().size();
    size_t n_vars = var_names.size();
    model.resize(n_cons, n_vars);

    Srm2Scm(model.getAr(), model.getAc());
}

void LpReader::parse_lincomb(const std::vector<std::string>& tokens, std::vector<Monom>& lhs,
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

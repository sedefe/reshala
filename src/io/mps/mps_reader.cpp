#include "reshala/io/mps/mps_reader.h"

#include <iostream>

namespace reshala {

FileReadStatus MpsReader::Read() {
    std::ifstream file(path_);
    if (!file.is_open()) {
        return FileReadStatus::kFsError;
    }

    std::string line;
    auto current_state = MpsParseState::kNon;
    model_.GetObj().orig_sense = Sense::kMin;

    while (std::getline(file, line)) {
        line_number++;
        while (line.back() == '\r') line.pop_back();

        auto tokens = tokenize_line(line);
        if (tokens.empty()) continue;

        if (tokens.size() == 2 and tokens[0] == "NAME") {
            current_state = MpsParseState::kNam;
            continue;
        } else if (tokens.size() == 1 and tokens[0] == "ROWS") {
            current_state = MpsParseState::kRow;
            continue;
        } else if (tokens.size() == 1 and tokens[0] == "COLUMNS") {
            // Now the size is known
            model_.GetAr().GetRows().resize(con_names.size(), 0);
            model_.GetRhs().resize(con_names.size());
            current_state = MpsParseState::kCol;
            continue;
        } else if (tokens.size() == 3 and tokens[2] == "'INTORG'") {
            int_marker = true;
            continue;
        } else if (tokens.size() == 3 and tokens[2] == "'INTEND'") {
            int_marker = false;
            continue;
        } else if (tokens.size() == 1 and tokens[0] == "RHS") {
            current_state = MpsParseState::kRhs;
            continue;
        } else if (tokens.size() == 1 and tokens[0] == "RANGES") {
            current_state = MpsParseState::kRhs;
            ThrowParseError("Don't support RANGES section in MPS now");
        } else if (tokens.size() == 1 and tokens[0] == "BOUNDS") {
            FinalizeRhs();
            size_t n_cons = con_names.size();
            size_t n_vars = var_names.size();
            model_.Resize(n_cons, n_vars);
            model_.FinalizeAc();

            current_state = MpsParseState::kBnd;
            continue;
        } else if (tokens.size() == 1 and tokens[0] == "ENDDATA") {
            current_state = MpsParseState::kDon;
            break;
        }

        switch (current_state) {
            case MpsParseState::kNam:
                break;
            case MpsParseState::kRow:
                ParseRows(tokens);
                break;
            case MpsParseState::kCol:
                ParseColumns(tokens);
                break;
            case MpsParseState::kRhs:
                ParseRhs(tokens);
                break;
            case MpsParseState::kBnd:
                ParseBounds(tokens);
                break;
            default:
                break;
        }
    }
    file.close();
    return FileReadStatus::kOk;
}

void MpsReader::ParseRows(const std::vector<std::string>& tokens) {
    if ((tokens.size() != 2) or (tokens[0].size() != 1)) {
        ThrowParseError("Bad ROW section");
    }
    char sense = tokens[0][0];
    const std::string& name = tokens[1];
    ExpType exp_type = MpsChar2ExpType(sense);

    if (exp_type != ExpType::kNon) {
        con_names.get_index(name);
        con_types.push_back(exp_type);
        con_rhs.push_back(0.0);
    } else {
        obj_name = name;
    }
}

void MpsReader::ParseColumns(const std::vector<std::string>& tokens) {
    if ((tokens.size() != 3) and (tokens.size() != 5)) {
        ThrowParseError("Bad COLUMNS section");
    }
    Index var_index = var_names.get_index(tokens[0]);
    if (var_index >= model_.GetVars().Size()) {
        model_.GetObj().coefficients.push_back(0.0);
        model_.GetVars().Push({}, int_marker);
    }
    for (Index i = 1; i < tokens.size(); i += 2) {
        Scalar coeff = std::stod(tokens[i + 1]);
        if (tokens[i] != obj_name) {
            Index con_index = con_names.get_index(tokens[i]);
            model_.GetAr().GetRow(con_index).Push(var_index, coeff);
        } else {
            model_.GetObj().coefficients[var_index] = coeff;
        }
    }
}

void MpsReader::ParseRhs(const std::vector<std::string>& tokens) {
    if (((tokens.size() != 3) and (tokens.size() != 5)) or (tokens[0] != "RHS")) {
        ThrowParseError("Bad RHS section");
    }
    for (Index i = 1; i < tokens.size(); i += 2) {
        Index con_index = con_names.get_index(tokens[i]);
        Scalar coeff = std::stod(tokens[i + 1]);
        con_rhs[con_index] = coeff;
    }
}

void MpsReader::FinalizeRhs() {
    for (Index ic = 0; ic < con_names.size(); ic++) {
        auto con_type = con_types[ic];
        auto rhs = con_rhs[ic];
        model_.GetRhs()[ic] = ExpType2Bounds(con_type, rhs);
    }
}

void MpsReader::ParseBounds(const std::vector<std::string>& tokens) {}

}  // namespace reshala

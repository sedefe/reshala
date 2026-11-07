#include "reshala/io/mps/mps_reader.h"

namespace reshala {

MpsBoundType Str2MpsBoundType(const std::string& s) {
    if (s == "LO") return MpsBoundType::kLO;
    if (s == "UP") return MpsBoundType::kUP;
    if (s == "FX") return MpsBoundType::kFX;
    if (s == "FR") return MpsBoundType::kFR;
    if (s == "MI") return MpsBoundType::kMI;
    if (s == "PL") return MpsBoundType::kPL;
    if (s == "BV") return MpsBoundType::kBV;
    if (s == "LI") return MpsBoundType::kLI;
    if (s == "UI") return MpsBoundType::kUI;
    return MpsBoundType::kNon;
}

FileReadStatus MpsReader::Parse(std::ifstream& file) {
    std::string line;
    auto current_state = MpsParseState::kNon;
    model_.GetObj().orig_sense = Sense::kMin;

    while (std::getline(file, line)) {
        line_number++;
        while (line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '*') continue;

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
            model_.Resize(GetNCons(), 0);
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
            current_state = MpsParseState::kRng;
            continue;
        } else if (tokens.size() == 1 and tokens[0] == "BOUNDS") {
            model_.Resize(GetNCons(), GetNVars());
            model_.FinalizeAc();

            current_state = MpsParseState::kBnd;
            continue;
        } else if (tokens.size() == 1 and tokens[0] == "ENDATA") {
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
            case MpsParseState::kRng:
                ParseRanges(tokens);
                break;
            case MpsParseState::kBnd:
                ParseBounds(tokens);
                break;
            default:
                ThrowParseError("Unexpected line: \"" + line + "\"");
        }
    }
    file.close();
    return (current_state == MpsParseState::kDon) ? FileReadStatus::kOk
                                                  : FileReadStatus::kParseError;
}

void MpsReader::ParseRows(const std::vector<std::string>& tokens) {
    if ((tokens.size() != 2) or (tokens[0].size() != 1)) {
        ThrowParseError("Bad ROW section");
    }
    char sense = tokens[0][0];
    const std::string& name = tokens[1];
    ExpType exp_type = MpsChar2ExpType(sense);

    if (exp_type != ExpType::kNon) {
        GetConIndex(name, true);
        con_types.push_back(exp_type);
        model_.GetRhs().push_back(ExpType2Bounds(exp_type, 0.0));
    } else {
        // "If more than one free row is specified, the first one is used as the objective function
        // and the others are discarded"
        if (!has_obj) {
            SetObjName(name);
            has_obj = true;
        } else {
            discarded_free_rows.insert(name);
        }
    }
}

void MpsReader::ParseColumns(const std::vector<std::string>& tokens) {
    if ((tokens.size() != 3) and (tokens.size() != 5)) {
        ThrowParseError("Bad COLUMNS section");
    }

    if (!GetNames().vars.Contains(tokens[0])) {  // Add a new one
        model_.GetObj().coefficients.push_back(0.0);
        if (int_marker) {
            // "If no bounds are specified for the variables within markers, bounds of 0 (zero) and
            // 1 (one) are assumed."
            model_.GetVars().Push({0., 1.}, true);
        } else {
            model_.GetVars().Push({}, false);
        }
    }
    Index var_index = GetVarIndex(tokens[0], true);
    for (Index i = 1; i < tokens.size(); i += 2) {
        if (discarded_free_rows.find(tokens[i]) != discarded_free_rows.end()) {
            continue;  // Skip free rows
        }

        Scalar coeff = std::stod(tokens[i + 1]);
        if (tokens[i] != GetObjName()) {
            Index ic = GetConIndex(tokens[i], false);
            if (ic >= model_.GetNCons()) {
                ThrowParseError("Constraint did not appear in ROWS: " + tokens[i]);
            }
            if (coeff != 0.0) {
                model_.GetRow(ic).Push(var_index, coeff);
            }
        } else {
            model_.GetObj().coefficients[var_index] = coeff;
        }
    }
}

void MpsReader::ParseRhs(const std::vector<std::string>& tokens) {
    if ((tokens.size() != 3) and (tokens.size() != 5)) {
        ThrowParseError("Bad RHS section");
    }
    for (Index i = 1; i < tokens.size(); i += 2) {
        if (tokens[i] == GetObjName()) {  // Todo: check specs for objective rhs
            continue;
        }
        Index ic = GetConIndex(tokens[i], false);
        Scalar rhs = std::stod(tokens[i + 1]);
        model_.GetRhs(ic) = ExpType2Bounds(con_types[ic], rhs);
    }
}

void MpsReader::ParseRanges(const std::vector<std::string>& tokens) {
    if ((tokens.size() != 3) and (tokens.size() != 5)) {
        ThrowParseError("Bad RANGES section");
    }
    for (Index i = 1; i < tokens.size(); i += 2) {
        Index ic = GetConIndex(tokens[i], false);
        Scalar range = std::stod(tokens[i + 1]);
        Bounds& rhs = model_.GetRhs(ic);
        switch (con_types[ic]) {
            case ExpType::kGe:
                model_.GetRhs(ic) = {rhs.le, rhs.le + std::abs(range)};
                break;
            case ExpType::kLe:
                model_.GetRhs(ic) = {rhs.ri - std::abs(range), rhs.ri};
                break;
            case ExpType::kEq:
                if (range > 0) {
                    model_.GetRhs(ic) = {rhs.le, rhs.ri + std::abs(range)};
                } else {
                    model_.GetRhs(ic) = {rhs.le - std::abs(range), rhs.ri};
                }
                break;
            default:
                ThrowParseError("Unknown type for constraint " + tokens[i]);
        }
    }
}

void MpsReader::ParseBounds(const std::vector<std::string>& tokens) {
    if ((tokens.size() != 3 and tokens.size() != 4) or tokens[0].size() != 2) {
        ThrowParseError("Bad BOUNDS section");
    }

    auto type = Str2MpsBoundType(tokens[0]);
    Index var_index = GetVarIndex(tokens[2], false);
    Scalar value;
    switch (type) {
        case MpsBoundType::kLI:  // integer kLO
            value = std::stod(tokens[3]);
            model_.SetIntegrality(var_index, true);
            if (value == 0.0) {
                // "To specify general integers with no upper bounds, use LI with the value 0.0."
                model_.SetBounds(var_index, {0, kInf});
            }
        case MpsBoundType::kLO:
            value = std::stod(tokens[3]);
            model_.SetBounds(var_index, {value, model_.GetBounds(var_index).ri});
            break;
        case MpsBoundType::kUI:  // integer kUI
            model_.SetIntegrality(var_index, true);
        case MpsBoundType::kUP:
            value = std::stod(tokens[3]);
            model_.SetBounds(var_index, {model_.GetBounds(var_index).le, value});
            break;
        case MpsBoundType::kFX:
            value = std::stod(tokens[3]);
            model_.SetBounds(var_index, {value, value});
            break;
        case MpsBoundType::kFR:
            model_.SetBounds(var_index, {-kInf, kInf});
            break;
        case MpsBoundType::kMI:
            model_.SetBounds(var_index, {-kInf, model_.GetBounds(var_index).ri});
            break;
        case MpsBoundType::kPL:
            model_.SetBounds(var_index, {model_.GetBounds(var_index).le, kInf});
            break;
        case MpsBoundType::kBV:
            model_.SetBounds(var_index, {0, 1});
            model_.SetIntegrality(var_index, true);
            break;
        default:
            ThrowParseError("Unsupported bound type: " + tokens[0]);
    }
}

}  // namespace reshala

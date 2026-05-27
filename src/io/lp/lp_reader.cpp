#include "reshala/io/lp/lp_reader.h"

namespace reshala {

// О том, что строка закончилась, мы точно сможем узнать только глядя на следующую строку. Поэтому
// для обжектива и ограничений суём все токены в вектор multiline и глядим на него только встретив
// знак сравнения (для ограничений) или встретив следующую секцию (для обжектива). В некоторых
// случаях секция ограничений опциональна, так что и эта логика не без огреха.

FileReadStatus LpReader::Read(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        return FileReadStatus::kFsError;
    }

    std::string line;
    std::string con_exp_token;
    auto current_state = LpParseState::kNon;

    while (std::getline(file, line)) {
        line_number++;
        while (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty() || line[0] == '\\') continue;

        auto tokens = tokenize_line(line);
        if (tokens.empty()) continue;
        std::string lower_line = to_lowercase(line);

        if (lower_line == "min" or lower_line == "minimize") {
            model_.GetObj().orig_sense = Sense::kMin;
            multiline = {};
            current_state = LpParseState::kObj;
            continue;
        } else if (lower_line == "max" or lower_line == "maximize") {
            model_.GetObj().orig_sense = Sense::kMax;
            multiline = {};
            current_state = LpParseState::kObj;
            continue;
        } else if (lower_line == "s.t." or lower_line == "subject to") {
            ParseObjective(multiline);
            multiline = {};
            current_state = LpParseState::kCon;
            continue;
        } else if (lower_line == "bounds") {
            if (!matrix_finalized) {
                FinalizeMatrix();
            }
            current_state = LpParseState::kBnd;
            continue;
        } else if (lower_line == "binaries") {
            if (!matrix_finalized) {
                FinalizeMatrix();
            }
            current_state = LpParseState::kBin;
            continue;
        } else if (lower_line == "generals" or lower_line == "integers") {
            if (!matrix_finalized) {
                FinalizeMatrix();
            }
            current_state = LpParseState::kGen;
            continue;
        } else if (lower_line == "end") {
            if (!matrix_finalized) {
                FinalizeMatrix();
            }
            current_state = LpParseState::kDon;
            break;
        }

        switch (current_state) {
            case LpParseState::kObj:
                multiline.insert(multiline.end(), tokens.begin(), tokens.end());
                break;
            case LpParseState::kCon:
                multiline.insert(multiline.end(), tokens.begin(), tokens.end());
                con_exp_token = tokens[tokens.size() - 2];
                if (con_exp_token[0] == '<' or con_exp_token[0] == '>' or con_exp_token[0] == '=') {
                    ParseConstraint(multiline);
                    multiline = {};
                }
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
    Index begin = 0;
    if (tokens.size() >= 1 and tokens[0].back() == ':') {
        names_.obj.assign(tokens[0], 0, tokens[0].length() - 1);
        begin = 1;
    }
    Index end = tokens.size();
    ParseLincomb(tokens, lhs, begin, end);

    Objective& obj = model_.GetObj();
    model_.GetObj().coefficients.resize(names_.vars.size());
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
    Index begin = 0;
    if (tokens.size() >= 1 and tokens[0].back() == ':') {
        std::string con_name(tokens[0].begin(), tokens[0].end() - 1);
        names_.cons.get_index(con_name);
        begin = 1;
    }
    Index end = tokens.size() - 2;

    ParseLincomb(tokens, lhs, begin, end);

    const std::string& exp_token = tokens[tokens.size() - 2];
    const std::string& rhs_token = tokens[tokens.size() - 1];

    if (!(exp_token[0] == '<' or exp_token[0] == '>' or exp_token[0] == '=') or
        !((exp_token.size() == 1) or (exp_token.size() == 2 and exp_token[1] == '='))) {
        ThrowParseError("Unsupported expression type: " + exp_token);
    }
    ExpType exp_type = LpChar2ExpType(exp_token[0]);
    Scalar coeff = std::stod(rhs_token);
    Bounds rhs = ExpType2Bounds(exp_type, coeff);

    SparseVector sv(names_.vars.size());
    sv.Reserve(lhs.size());
    for (const auto& m : lhs) {
        sv.Push(m.index, m.coeff);
    }
    model_.PrepareConstraint(sv, rhs);
}

void LpReader::ParseBounds(const std::vector<std::string>& tokens) {
    size_t n = tokens.size();
    if (n != 3 and n != 5) ThrowParseError("Can't parse bounds");

    Index i_var = n - 3;  // 0 for "x <= 1", 2 for "0 <= x <= 1"
    if (names_.vars.name_to_index.find(tokens[i_var]) == names_.vars.name_to_index.end()) {
        ThrowParseError("Unexpected variable " + tokens[i_var] + " in Bounds section");
    }
    Index index = names_.vars.get_index(tokens[i_var]);

    std::vector<std::tuple<bool, ExpType, Scalar>> expressions;
    if (n == 3) {
        expressions.push_back({false, LpChar2ExpType(tokens[1][0]), std::stod(tokens[2])});
    } else {  // (n == 5)
        expressions.push_back({true, LpChar2ExpType(tokens[1][0]), std::stod(tokens[0])});
        expressions.push_back({false, LpChar2ExpType(tokens[3][0]), std::stod(tokens[4])});
    }

    for (auto& tuple : expressions) {
        bool do_swap = std::get<0>(tuple);
        ExpType type = std::get<1>(tuple);
        Scalar rhs = std::get<2>(tuple);
        if (do_swap) {  // "0 <= x" -> "x >= 0"
            if (type == ExpType::kGe)
                type = ExpType::kLe;
            else if (type == ExpType::kLe)
                type = ExpType::kGe;
        }
        const Bounds& bnd = model_.GetBounds(index);
        switch (type) {
            case ExpType::kGe:
                model_.SetBounds(index, BoundsIntersection(bnd, {rhs, kInf}));
                break;
            case ExpType::kLe:
                model_.SetBounds(index, BoundsIntersection(bnd, {-kInf, rhs}));
                break;
            case ExpType::kEq:
                model_.SetBounds(index, BoundsIntersection(bnd, {rhs, rhs}));
                break;
            default:
                assert(false && "Unknown ExpType");
        }
    }
}

void LpReader::ParseBinaries(const std::vector<std::string>& tokens) {
    for (const auto& str : tokens) {
        if (names_.vars.name_to_index.find(str) == names_.vars.name_to_index.end()) {
            ThrowParseError("Unexpected variable " + str + " in Binaries section");
        }
        auto index = names_.vars.get_index(str);
        model_.SetIntegrality(index, true);
        const Bounds& bnd = model_.GetBounds(index);
        model_.SetBounds(index, BoundsIntersection(bnd, {0, 1}));
    }
}

void LpReader::ParseGenerals(const std::vector<std::string>& tokens) {
    for (const auto& str : tokens) {
        if (names_.vars.name_to_index.find(str) == names_.vars.name_to_index.end()) {
            ThrowParseError("Unexpected variable " + str + " in Generals section");
        }
        auto index = names_.vars.get_index(str);
        model_.SetIntegrality(index, true);
    }
}

void sort_monoms(std::vector<Monom>& monoms) {
    std::sort(monoms.begin(), monoms.end(),
              [](const Monom& a, const Monom& b) { return a.index < b.index; });
}

void LpReader::ParseLincomb(const std::vector<std::string>& tokens, std::vector<Monom>& lhs,
                            Index begin, Index end) {
    Scalar sign = 1.0;
    Scalar coeff = 1.0;

    Index i;
    for (i = begin; i < end; ++i) {
        const std::string& token = tokens[i];

        if (token == "+") {
            sign = 1.0;
        } else if (token == "-") {
            sign = -1.0;
        } else if (isdigit(token[0]) || token[0] == '-' || token[0] == '.') {
            coeff = std::stod(token);
        } else {
            size_t idx = names_.vars.get_index(token);
            lhs.push_back({sign * coeff, idx});

            sign = 1.0;
            coeff = 1.0;
        }
    }
    sort_monoms(lhs);
}

void LpReader::FinalizeMatrix() {
    size_t n_cons = model_.GetAr().GetRows().size();
    size_t n_vars = names_.vars.size();
    model_.Resize(n_cons, n_vars);
    model_.FinalizeAc();
    matrix_finalized = true;
}

}  // namespace reshala

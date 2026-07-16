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
    Scalar c0;
    Index begin = 0;
    if (tokens.size() >= 1 and tokens[0].back() == ':') {
        names_.obj.assign(tokens[0], 0, tokens[0].length() - 1);
        begin = 1;
    }
    Index end = tokens.size();
    ParseLincomb(tokens, lhs, c0, begin, end);

    Objective& obj = model_.GetObj();
    model_.GetObj().coefficients.resize(names_.vars.Size());
    obj.c0 = c0;
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
    Scalar c0;
    Index begin = 0;
    if (tokens.size() >= 1 and tokens[0].back() == ':') {
        std::string con_name(tokens[0].begin(), tokens[0].end() - 1);
        names_.cons.GetIndex(con_name);
        begin = 1;
    }
    Index end = tokens.size() - 2;

    ParseLincomb(tokens, lhs, c0, begin, end);
    assert(c0 == 0.0 && "Free terms in constraints are not supported");

    const std::string& exp_token = tokens[tokens.size() - 2];
    const std::string& rhs_token = tokens[tokens.size() - 1];

    if (!(exp_token[0] == '<' or exp_token[0] == '>' or exp_token[0] == '=') or
        !((exp_token.size() == 1) or (exp_token.size() == 2 and exp_token[1] == '='))) {
        ThrowParseError("Unsupported expression type: " + exp_token);
    }
    ExpType exp_type = LpChar2ExpType(exp_token[0]);
    Scalar coeff = std::stod(rhs_token);
    Bounds rhs = ExpType2Bounds(exp_type, coeff);

    SparseVector sv(names_.vars.Size());
    sv.Reserve(lhs.size());
    for (const auto& m : lhs) {
        sv.Push(m.index, m.coeff);
    }
    model_.PrepareConstraint(sv, rhs);
}

void LpReader::ParseBounds(const std::vector<std::string>& tokens) {
    size_t n = tokens.size();
    if (n == 2 and tokens[1] == "Free") {
        Index index = names_.vars.GetIndex(tokens[0]);
        model_.SetBounds(index, {-kInf, kInf});
        return;
    }
    if (n != 3 and n != 5) ThrowParseError("Can't parse bounds");

    Index index;
    std::vector<std::tuple<bool, ExpType, Scalar>> expressions;
    if (n == 3) {
        bool var_at_0 = names_.vars.Contains(tokens[0]);
        bool var_at_2 = names_.vars.Contains(tokens[2]);
        if ((var_at_0 ^ var_at_2) == false) ThrowParseError("Can't parse bounds");

        if (var_at_0) {  // "x >= 0"  ->  no swap
            index = names_.vars.GetIndex(tokens[0]);
            expressions.push_back({false, LpChar2ExpType(tokens[1][0]), std::stod(tokens[2])});
        } else {  // "0 <= x"  ->  swap
            index = names_.vars.GetIndex(tokens[2]);
            expressions.push_back({true, LpChar2ExpType(tokens[1][0]), std::stod(tokens[0])});
        }
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
        auto index = names_.vars.GetIndex(str);
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
        auto index = names_.vars.GetIndex(str);
        model_.SetIntegrality(index, true);
    }
}

void sort_monoms(std::vector<Monom>& monoms) {
    std::sort(monoms.begin(), monoms.end(),
              [](const Monom& a, const Monom& b) { return a.index < b.index; });
}

// <expression> ::= [sign] <monom> { ('+' | '-') <monom> } [ ('+' | '-') <constant> ]
// <monom>      ::= [number] <variable>
// <constant>   ::= [number]   (only allowed as the very last term)
void LpReader::ParseLincomb(const std::vector<std::string>& tokens, std::vector<Monom>& lhs,
                            Scalar& free_term, Index begin, Index end) {
    enum State { kExpectOperand, kExpectVariable, kExpectOperator };
    State state = kExpectOperand;
    Scalar sign = 1.0;        // sign of the current monom
    Scalar coeff = 1.0;       // coefficient (if a number was seen)
    bool have_sign = false;   // whether we already set a sign for this monom
    bool has_number = false;  // whether we already parsed a number

    if (begin == end) return;  // Empty one
    for (Index i = begin; i < end; ++i) {
        const std::string& token = tokens[i];

        if (token == "+" || token == "-") {
            if (state == kExpectOperator) {
                // Valid operator between monoms – start a new monom
                sign = (token == "+") ? 1.0 : -1.0;
                have_sign = true;
                coeff = 1.0;
                has_number = false;
                state = kExpectOperand;
            } else if (state == kExpectOperand) {
                // Sign at the beginning of a monom (allowed only once)
                if (have_sign) ThrowParseError("Double operator at token " + std::to_string(i));
                sign = (token == "+") ? 1.0 : -1.0;
                have_sign = true;
                // remain in kExpectOperand, still need a variable or number
            } else {
                ThrowParseError("Token " + std::to_string(i) + ": expected variable name, got \"" +
                                token + "\"");
            }
        } else if (isdigit(token[0]) || token[0] == '-' || token[0] == '.') {
            // Number token
            // Allow implicit '+' when a number follows a variable
            if (state == kExpectOperator) {
                // Implicit plus: start a new monom
                sign = 1.0;
                have_sign = false;
                coeff = 1.0;
                has_number = false;
                state = kExpectOperand;
            }
            if (state != kExpectOperand) {
                ThrowParseError("Token " + std::to_string(i) + ": expected operand, got \"" +
                                token + "\"");
            }
            coeff = std::stod(token);
            has_number = true;
            state = kExpectVariable;
        } else {
            // Variable token
            if (state == kExpectOperator)
                ThrowParseError("Token " + std::to_string(i) + ": expected operator, got \"" +
                                token + "\"");

            Index idx = names_.vars.GetIndex(token);
            lhs.push_back({sign * coeff, idx});

            // Reset
            sign = 1.0;
            coeff = 1.0;
            have_sign = false;
            has_number = false;
            state = kExpectOperator;
        }
    }

    free_term = 0.0;
    if (state == kExpectVariable && has_number) {
        free_term = sign * coeff;
        state = kExpectOperator;
    }

    if (state != kExpectOperator) ThrowParseError("Unexpected end of a linear combination");

    sort_monoms(lhs);
}

void LpReader::FinalizeMatrix() {
    size_t n_cons = model_.GetAr().GetRows().size();
    size_t n_vars = names_.vars.Size();
    model_.Resize(n_cons, n_vars);
    model_.FinalizeAc();
    matrix_finalized = true;
}

}  // namespace reshala

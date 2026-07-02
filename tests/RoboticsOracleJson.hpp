#pragma once

// ============================================================================
//  RoboticsOracleJson.hpp -- the minimal JSON DOM + recursive-descent parser
//  shared by every robotics-oracle fixture loader (Scope A:
//  RoboticsOracleLoader.hpp; Scope B: RoboticsOracleMoleculeLoader.hpp). Split
//  out of RoboticsOracleLoader.hpp (Rule 2/3: Scope B's ".moldyn.manifest.json"
//  reader needs the exact same tiny parser -- duplicating ~150 lines across
//  two loader headers is worse than one shared header; behavior is otherwise
//  byte-for-byte what TestRoboticsOracle.cpp already exercised, so this is a
//  pure extraction, not a behavior change).
//
//  JSON: every manifest here is small, fixed-schema, machine-written-by-us-
//  only JSON (docs/specs/robotics-oracle-differential.md §9: "small and
//  human-reviewable"). A full JSON library dependency buys nothing this
//  ~150-line parser doesn't already give (Rule 2) -- it handles the JSON
//  subset our generators emit (objects/arrays/strings/numbers/bool/null),
//  nothing more.
// ============================================================================

#include <cctype>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace robotics_oracle_loader {

// ---------------------------------------------------------------------------
//  Minimal JSON DOM + recursive-descent parser (see banner).
// ---------------------------------------------------------------------------
struct JsonValue {
    enum class Type { Null, Bool, Number, String, Array, Object } type = Type::Null;
    bool boolVal = false;
    double numVal = 0;
    std::string strVal;
    std::vector<JsonValue> arrVal;
    std::vector<std::pair<std::string, JsonValue>> objVal;

    const JsonValue& at(const std::string& key) const {
        for (const auto& kv : objVal) {
            if (kv.first == key) {
                return kv.second;
            }
        }
        throw std::runtime_error("robotics_oracle: manifest JSON missing key '" + key + "'");
    }
    bool has(const std::string& key) const {
        for (const auto& kv : objVal) {
            if (kv.first == key) {
                return true;
            }
        }
        return false;
    }
};

class JsonParser {
public:
    explicit JsonParser(const std::string& s) : s_(s) {}
    JsonValue parse() {
        skipWs();
        return parseValue();
    }

private:
    const std::string& s_;
    std::size_t i_ = 0;

    void skipWs() {
        while (i_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[i_]))) {
            ++i_;
        }
    }
    char peek() const { return i_ < s_.size() ? s_[i_] : '\0'; }

    JsonValue parseValue() {
        skipWs();
        const char c = peek();
        if (c == '{') return parseObject();
        if (c == '[') return parseArray();
        if (c == '"') return parseString();
        if (c == 't' || c == 'f') return parseBool();
        if (c == 'n') return parseNull();
        return parseNumber();
    }
    JsonValue parseObject() {
        JsonValue v;
        v.type = JsonValue::Type::Object;
        ++i_;
        skipWs();
        if (peek() == '}') {
            ++i_;
            return v;
        }
        while (true) {
            skipWs();
            JsonValue key = parseString();
            skipWs();
            if (peek() != ':') {
                throw std::runtime_error("robotics_oracle: manifest JSON: expected ':'");
            }
            ++i_;
            JsonValue val = parseValue();
            v.objVal.emplace_back(key.strVal, std::move(val));
            skipWs();
            if (peek() == ',') {
                ++i_;
                continue;
            }
            if (peek() == '}') {
                ++i_;
                break;
            }
            throw std::runtime_error("robotics_oracle: manifest JSON: expected ',' or '}'");
        }
        return v;
    }
    JsonValue parseArray() {
        JsonValue v;
        v.type = JsonValue::Type::Array;
        ++i_;
        skipWs();
        if (peek() == ']') {
            ++i_;
            return v;
        }
        while (true) {
            v.arrVal.push_back(parseValue());
            skipWs();
            if (peek() == ',') {
                ++i_;
                continue;
            }
            if (peek() == ']') {
                ++i_;
                break;
            }
            throw std::runtime_error("robotics_oracle: manifest JSON: expected ',' or ']'");
        }
        return v;
    }
    JsonValue parseString() {
        if (peek() != '"') {
            throw std::runtime_error("robotics_oracle: manifest JSON: expected string");
        }
        ++i_;
        JsonValue v;
        v.type = JsonValue::Type::String;
        std::string out;
        while (i_ < s_.size() && s_[i_] != '"') {
            char c = s_[i_];
            if (c == '\\' && i_ + 1 < s_.size()) {
                ++i_;
                const char e = s_[i_];
                switch (e) {
                    case 'n': out.push_back('\n'); break;
                    case 't': out.push_back('\t'); break;
                    case '"': out.push_back('"'); break;
                    case '\\': out.push_back('\\'); break;
                    case '/': out.push_back('/'); break;
                    default: out.push_back(e); break;
                }
            } else {
                out.push_back(c);
            }
            ++i_;
        }
        if (i_ >= s_.size()) {
            throw std::runtime_error("robotics_oracle: manifest JSON: unterminated string");
        }
        ++i_; // closing quote
        v.strVal = out;
        return v;
    }
    JsonValue parseBool() {
        JsonValue v;
        v.type = JsonValue::Type::Bool;
        if (s_.compare(i_, 4, "true") == 0) {
            v.boolVal = true;
            i_ += 4;
        } else if (s_.compare(i_, 5, "false") == 0) {
            v.boolVal = false;
            i_ += 5;
        } else {
            throw std::runtime_error("robotics_oracle: manifest JSON: bad literal");
        }
        return v;
    }
    JsonValue parseNull() {
        JsonValue v;
        v.type = JsonValue::Type::Null;
        if (s_.compare(i_, 4, "null") == 0) {
            i_ += 4;
        } else {
            throw std::runtime_error("robotics_oracle: manifest JSON: bad literal");
        }
        return v;
    }
    JsonValue parseNumber() {
        const std::size_t start = i_;
        if (peek() == '-') {
            ++i_;
        }
        while (i_ < s_.size() &&
               (std::isdigit(static_cast<unsigned char>(s_[i_])) || s_[i_] == '.' || s_[i_] == 'e' ||
                s_[i_] == 'E' || s_[i_] == '+' || s_[i_] == '-')) {
            ++i_;
        }
        JsonValue v;
        v.type = JsonValue::Type::Number;
        v.numVal = std::stod(s_.substr(start, i_ - start));
        return v;
    }
};

inline std::string readWholeTextFile(const std::string& path) {
    std::ifstream f(path);
    if (!f) {
        throw std::runtime_error("robotics_oracle: cannot open '" + path + "'");
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

inline JsonValue loadManifestJsonFile(const std::string& path) {
    const std::string text = readWholeTextFile(path);
    JsonParser p(text);
    return p.parse();
}

inline JsonValue loadManifestJson(const std::string& fixturesDir, const std::string& caseName) {
    return loadManifestJsonFile(fixturesDir + "/" + caseName + ".manifest.json");
}

// Process-lifetime interning for the schema's `const char*` label/name
// fields (RoboticsOracleTypes.hpp predates this loader and is engine-
// agnostic/unchanged, Rule 3 -- it still expects string-literal-lifetime
// pointers). Loaded fixtures live for the whole test binary's run, so a
// simple never-freed pool is sufficient and safe (single-threaded loading).
inline const char* internString(const std::string& s) {
    static std::vector<std::unique_ptr<std::string>> pool;
    pool.push_back(std::make_unique<std::string>(s));
    return pool.back()->c_str();
}

} // namespace robotics_oracle_loader

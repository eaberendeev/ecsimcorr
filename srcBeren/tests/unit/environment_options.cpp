#include "env_options.h"
#include "test_support.h"

void environment_options() {
    const char* name = "BEREN_TEST_ENV_PARSER";
    check::equal(envOptions::detail::parseIntegerEnv<int>(name, "+42"), 42, "integer plus sign");
    check::equal(envOptions::detail::parseIntegerEnv<int>(name, "-42"), -42, "integer minus sign");
    check::equal(envOptions::detail::parseIntegerEnv<int64_t>(name, "9223372036854775807"),
                 std::numeric_limits<int64_t>::max(), "int64 maximum");
    for (const char* value : {"", " ", "+", "+-1", "++1", "1x", " 1", "1 ", "1\t", "2147483648", "-2147483649"}) {
        bool rejected = false;
        try {
            envOptions::detail::parseIntegerEnv<int>(name, value);
        } catch (const std::invalid_argument& error) {
            rejected = std::string(error.what()).find(name) != std::string::npos;
        }
        check::require(rejected, "invalid integer must identify its environment variable");
    }
    const char* previous = std::getenv(name);
    const bool was_set = previous != nullptr;
    const std::string saved = was_set ? previous : "";
    for (const char* value : {"true", "YES", "1", "false", "No", "0"}) {
        setenv(name, value, 1);
        const bool expected = std::string(value) == "true" || std::string(value) == "YES" || std::string(value) == "1";
        check::equal(envOptions::getenvParsed<bool>(name, !expected), expected, "boolean parsing");
    }
    for (const char* value : {"", " true", "false ", "1\t", "maybe"}) {
        setenv(name, value, 1);
        bool rejected = false;
        try {
            envOptions::getenvParsed<bool>(name, false);
        } catch (const std::invalid_argument& error) {
            rejected = std::string(error.what()).find(name) != std::string::npos;
        }
        check::require(rejected, "invalid boolean must identify its environment variable");
    }
    unsetenv(name);
    check::require(envOptions::getenvParsed<bool>(name, true), "unset boolean uses default");
    check::equal(envOptions::getenvParsed<int>(name, 17), 17, "unset integer uses default");
    if (was_set)
        setenv(name, saved.c_str(), 1);
}

int main() {
    return check::run(environment_options);
}

#pragma once

#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

template <typename T>
T get_checked(const nlohmann::json& j, const std::string& key) {
    if (!j.contains(key)) {
        throw std::runtime_error("Missing key: " + key);
    }
    return j[key].get<T>();
}

// Специализации для более информативных сообщений
template <>
inline std::string get_checked<std::string>(const nlohmann::json& j, const std::string& key) {
    if (!j.contains(key)) {
        throw std::runtime_error("Missing key: " + key);
    }

    if (!j[key].is_string()) {
        throw std::runtime_error("Key '" + key + "' is not a string");
    }

    return j[key].get<std::string>();
}

template <>
inline int get_checked<int>(const nlohmann::json& j, const std::string& key) {
    if (!j.contains(key)) {
        throw std::runtime_error("Missing key: " + key);
    }

    if (!j[key].is_number_integer()) {
        throw std::runtime_error("Key '" + key + "' is not an integer");
    }

    return j[key].get<int>();
}

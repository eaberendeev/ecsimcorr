// // Author: Evgeny Berendeev
// // Email: evgeny.berendeev@gmail.com
// // Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>

inline constexpr int double_to_int(const double d) {
    return static_cast<int>(d + 1.0) - 1;
}

inline int ngp(const double normalized_coord) {
    return std::lround(normalized_coord);
}

bool create_directory(const std::string& path);
std::vector<std::string> split_string(const std::string& s, const char delim);

inline int int_value(const double dblValue) {
    return static_cast<int>(dblValue < 0 ? dblValue - 0.5 : dblValue + 0.5);
}

template <typename T>
inline int sign(const T val) {
    return (T(0) < val) - (val < T(0));
}

inline std::string to_string(const int i, const int len) {
    std::stringstream ss;
    ss << std::setw(len) << std::setfill('0') << i;
    std::string s = ss.str();
    return s;
}

// inline Vector3i cast_to_int(const Vector3d& value) {
//     return Vector3i(static_cast<int>(value(0)), static_cast<int>(value(1)),
//                     static_cast<int>(value(2)));
// }

// inline Vector3i vector_to_vector3i(const std::vector<int>& vec) {
//     assert(vec.size() == 3 &&
//            "Std::vector has invalid size. Convertation to Vecto3d: FAILED");
//     return Vector3i(vec[0], vec[1], vec[2]);
// }
// inline Vector3d vector_to_vector3d(const std::vector<double>& vec) {
//     assert(vec.size() == 3 &&
//            "Std::vector has invalid size. Convertation to Vecto3d: FAILED");
//     return Vector3d(vec[0], vec[1], vec[2]);
// }

inline std::vector<double> string_to_doubles(std::vector<std::string> words) {
    std::vector<double> doubles;
    for (const auto& word : words) {
        doubles.push_back(std::stod(word));
    }
    return doubles;
}

template <typename Func>
std::pair<double, double> find_maximum_universal(Func f, double x_min, double x_max, int intervals) {
    // if (intervals <= 0 || x_min >= x_max) {
    //     throw std::invalid_argument("Invalid input parameters");
    // }

    double step = (x_max - x_min) / intervals;
    double max_value = f(x_min);
    double x_optimal = x_min;

    for (int i = 1; i <= intervals; ++i) {
        double x = x_min + i * step;
        double current_value = f(x);
        if (current_value > max_value) {
            max_value = current_value;
            x_optimal = x;
        }
    }

    return {max_value, x_optimal};
}

inline int hash(const std::string& key, int tableSize) {
    int hashVal = 0;

    for (size_t i = 0; i < key.length(); i++) {
        hashVal = 37 * hashVal + key[i];
    }

    hashVal %= tableSize;

    if (hashVal < 0)
        hashVal += tableSize;

    return hashVal;
}

struct pair_hash {
    size_t operator()(const std::pair<int, int>& p) const {
        return static_cast<size_t>(p.first) << 32 | p.second;
    }
};

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

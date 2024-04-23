// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef PARAMETERS_MAP_H
#define PARAMETERS_MAP_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

class ParametersMap {
   public:
    void set(const std::string& key, const std::vector<std::string>& values) {
        vectorMap.insert({key, values});
    }

    double get_double(const std::string& key, const int pos = 0) const {
        if (vectorMap.count(key) == 0) {
            std::cout << "key " << key << " not found\n";
        }
        return stod(get_values(key).at(pos));
    }

    int get_int(const std::string& key, const int pos = 0) const {
        if (vectorMap.count(key) == 0) {
           std::cout << "key " << key << " not found\n";
        }
        return stoi(get_values(key).at(pos));
    }

    const std::vector<std::string>& get_values(const std::string& key) const {
        return vectorMap.at(key);
    }

    bool is_empty() const { return vectorMap.empty(); }

    void print() const{
        for(const auto& elem: vectorMap){
            std::cout << elem.first << ": {";
            for(const auto& value : elem.second){
                std::cout << value << ", ";
            } 
            std::cout << "}\n";
        }
    }

   private:
    std::map<std::string, std::vector<std::string>> vectorMap;
};

ParametersMap load_parameters(const std::string& nameOfFileParameters);
std::vector<ParametersMap> load_vector_parameters(
    const std::string& nameOfFileParameters, const std::string& parametersDelim);

#endif

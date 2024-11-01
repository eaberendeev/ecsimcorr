// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "parameters_map.h"

#include <string.h>

#include <fstream>
#include <iostream>

#include "service.h"

ParametersMap load_parameters(const std::string& nameOfFileParameters) {
    ParametersMap parametersMap;
    std::vector<std::string> strvec;
    std::ifstream in(nameOfFileParameters);
    std::string line;

    std::string bufKey, bufVal;

    while (in.good()) {
        getline(in, line);

        if (in.eof())
            break;
        // only parameters may be set
        if (line.size() == 0)
            continue;

        strvec = split_string(line, ' ');

        std::vector<std::string>::iterator it = remove_if(
            strvec.begin(), strvec.end(), std::mem_fn(&std::string::empty));
        strvec.erase(it, strvec.end());

        if (strvec.size() < 2)
            continue;
        if (strvec[0][0] == '#')
            continue;
        // get name of parameter
        bufKey = strvec[0];
        // get values of parameters
        std::vector<std::string> bufVals;
        for (std::size_t i = 1; i < strvec.size(); i++) {
            bufVals.push_back(strvec[i]);
        }
        parametersMap.set(bufKey, bufVals);
    }
    return parametersMap;
}

std::vector<ParametersMap> load_vector_parameters(
    const std::string& nameOfFileParameters,
    const std::string& parametersDelim) {
    std::vector<ParametersMap> vectorParameters;
    std::vector<std::string> strvec;
    std::ifstream in(nameOfFileParameters);
    std::string line;

    std::string bufKey, bufVal;
    int currentParametr = -1;

    while (in.good()) {
        getline(in, line);

        if (in.eof())
            break;
        // only parameters may be set
        if (line.size() == 0)
            continue;

        strvec = split_string(line, ' ');

        std::vector<std::string>::iterator it = remove_if(
            strvec.begin(), strvec.end(), std::mem_fn(&std::string::empty));
        strvec.erase(it, strvec.end());

        if (strvec.size() < 2)
            continue;
        if (strvec[0][0] == '#')
            continue;
        // get name of parameter
        bufKey = strvec[0];
        if (bufKey == parametersDelim) {
            ParametersMap parametrsMap;
            vectorParameters.push_back(parametrsMap);
            currentParametr++;
        }
        if(currentParametr < 0) {
            std::cerr<< "Warning! Load vector parameters from file " << nameOfFileParameters << 
            " can be incorrect. Delimiter '" << parametersDelim << "' not found in header!\n"; 
        continue;
    }
        // get values of parameters
        std::vector<std::string> bufVals;
        for (std::size_t i = 1; i < strvec.size(); i++) {
            bufVals.push_back(strvec[i]);
        }
        vectorParameters.at(currentParametr).set(bufKey, bufVals);
    }
    return vectorParameters;
}

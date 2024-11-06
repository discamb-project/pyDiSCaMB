#pragma once

#include <string>
#include <map>
#include <vector>

struct GaussianScatteringParameters {
    std::vector<double> a;
    std::vector<double> b;
    double c;
    std::string repr();
};

std::map<std::string, GaussianScatteringParameters> get_table(std::string table);

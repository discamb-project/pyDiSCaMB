#pragma once

#include <string>
#include <map>
#include <vector>

struct GaussianScatterer {
    std::vector<double> a;
    std::vector<double> b;
    double c;
    std::string repr();
};

std::map<std::string, GaussianScatterer> get_table(std::string table);

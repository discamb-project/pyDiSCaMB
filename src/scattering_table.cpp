#include "scattering_table.hpp"

#include "discamb/Scattering/NGaussianFormFactorsTable.h"
#include "discamb/BasicChemistry/periodic_table.h"

#include <sstream>
#include <iomanip>

using namespace std;
using namespace discamb;

string GaussianScatterer::repr(){
    stringstream out;
    out << fixed << setprecision(2);
    out << "GaussianScatterer(" << endl;
    out << "    a = [";
    for (double ai: a){
        out << ai << ", ";
    }
    out << "]," << endl;
    out << "    b = [";
    for (double bi: b){
        out << bi << ", ";
    }
    out << "]," << endl;
    out << "    c = " << c << endl;
    out << ")";
    return out.str();
}

vector<string> get_possible_entries(string atom){
    vector<string> out {atom};
    int charge;
    // Positive
    for (charge = 1; charge <= 6; charge++){
        out.push_back(atom + to_string(charge) + string("+"));
    }
    // Negative
    for (charge = 1; charge <= 2; charge++){
        out.push_back(atom + to_string(charge) + string("-"));
    }
    out.push_back(atom + string("val"));
    out.push_back(atom + string("iso"));
    return out;
}

map<string, GaussianScatterer> get_table(string table){
    map<string, GaussianScatterer> out;
    for (int z = 1; z < 120; z++){
        string atom = periodic_table::symbol(z);
        vector<string> possible_entries = get_possible_entries(atom);
        for (string possible_entry : possible_entries){
            if (!n_gaussian_form_factors_table::hasFormFactor(possible_entry, table)){
                continue;
            }
            NGaussianFormFactor ff = n_gaussian_form_factors_table::getFormFactor(possible_entry, table);
            GaussianScatterer s;
            ff.get_parameters(s.a, s.b, s.c);
            out[possible_entry] = s;
        }
    }
    return out;
}

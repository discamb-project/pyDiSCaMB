#include "scattering_table.hpp"

#include <iomanip>
#include <sstream>

#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/Scattering/NGaussianFormFactorsTable.h"

using namespace std;
using namespace discamb;

string table_alias(string table) {
    if (table == "it") return "IT92";
    if (table == "IT") return "IT92";
    if (table == "it92") return "IT92";
    if (table == "it1992") return "IT92";
    if (table == "IT1992") return "IT92";
    if (table == "it-92") return "IT92";
    if (table == "it-1992") return "IT92";
    if (table == "IT-92") return "IT92";
    if (table == "IT-1992") return "IT92";

    if (table == "wk") return "Waasmeier-Kirfel";
    if (table == "WK") return "Waasmeier-Kirfel";
    if (table == "wk95") return "Waasmeier-Kirfel";
    if (table == "wk1995") return "Waasmeier-Kirfel";
    if (table == "WK95") return "Waasmeier-Kirfel";
    if (table == "WK1995") return "Waasmeier-Kirfel";
    if (table == "wk-95") return "Waasmeier-Kirfel";
    if (table == "wk-1995") return "Waasmeier-Kirfel";
    if (table == "WK-95") return "Waasmeier-Kirfel";
    if (table == "WK-1995") return "Waasmeier-Kirfel";
    if (table == "Waasmeier-Kirfel95") return "Waasmeier-Kirfel";
    if (table == "Waasmeier-Kirfel1995") return "Waasmeier-Kirfel";
    if (table == "Waasmeier-Kirfel-95") return "Waasmeier-Kirfel";
    if (table == "Waasmeier-Kirfel-1995") return "Waasmeier-Kirfel";

    if (table == "xray") return "IT92";

    if (table == "electron") return "electron-cctbx";

    if (table == "n_gaussian") return "IT92";

    return table;
}

vector<string> get_possible_entries(string atom) {
    vector<string> out{atom};
    int charge;
    // Positive
    for (charge = 1; charge <= 6; charge++) {
        out.push_back(atom + to_string(charge) + string("+"));
    }
    // Negative
    for (charge = 1; charge <= 2; charge++) {
        out.push_back(atom + to_string(charge) + string("-"));
    }
    out.push_back(atom + string("val"));
    out.push_back(atom + string("iso"));
    return out;
}

map<string, GaussianScatteringParameters> get_table(string table) {
    string alias = table_alias(table);
    map<string, GaussianScatteringParameters> out;
    for (int z = 1; z < 120; z++) {
        string atom = periodic_table::symbol(z);
        vector<string> possible_entries = get_possible_entries(atom);
        for (string possible_entry : possible_entries) {
            if (!n_gaussian_form_factors_table::hasFormFactor(possible_entry,
                                                              alias)) {
                continue;
            }
            NGaussianFormFactor ff =
                n_gaussian_form_factors_table::getFormFactor(possible_entry,
                                                             alias);
            GaussianScatteringParameters s;
            ff.get_parameters(s.a, s.b, s.c);
            out[possible_entry] = s;
        }
    }
    return out;
}

string GaussianScatteringParameters::repr() {
    stringstream out;
    out << fixed << setprecision(2);
    out << "GaussianScatteringParameters(" << endl;
    out << "    a = [";
    for (double ai : a) {
        out << ai << ", ";
    }
    out << "]," << endl;
    out << "    b = [";
    for (double bi : b) {
        out << bi << ", ";
    }
    out << "]," << endl;
    out << "    c = " << c << endl;
    out << ")";
    return out.str();
}

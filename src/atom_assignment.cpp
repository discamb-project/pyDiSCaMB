// This file is a modification of DiSCaMB's typeQuest example.
// main has been renamed to write_assignment_logs, and
// readStructureAndAssign has been renamed to read_bank_and_assign_atoms,
// both with different input parameters (instead of creating them inside the
// funcitons). File input features are also removed.

#include "atom_assignment.hpp"

#include <chrono>
#include <ctime>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/AtomTyping/atom_typing_utilities.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/discamb_version.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/HC_Model/hc_model_utilities.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/MathUtilities/graph_algorithms.h"
#include "discamb/MathUtilities/real_spherical_harmonics.h"
#include "discamb/MathUtilities/statistics.h"
#include "discamb/StructuralProperties/structural_properties.h"

using namespace discamb;
using namespace std;

string atomicNumbers2Formula(const vector<int> &atomicNumbers) {
    map<string, int> formula;
    string symbol, formulaStr;
    for (auto z : atomicNumbers) {
        symbol = periodic_table::symbol(z);
        if (formula.find(symbol) != formula.end())
            formula[symbol]++;
        else
            formula[symbol] = 1;
    }
    for (auto x : formula) {
        formulaStr += x.first;
        if (x.second != 1) formulaStr += to_string(x.second);
    }
    return formulaStr;
}

struct Atoms1stNeighbours {
    string symbol;
    string formula;
    bool operator>(const Atoms1stNeighbours &x) {
        pair<string, string> p1(symbol, formula), p2(x.symbol, x.formula);
        return p1 > p2;
    }
};

bool operator<(const Atoms1stNeighbours &x1, const Atoms1stNeighbours &x2) {
    pair<string, string> p1(x1.symbol, x1.formula), p2(x2.symbol, x2.formula);
    return p1 < p2;
}

bool less1(const Atoms1stNeighbours &v1, const Atoms1stNeighbours &v2) {
    pair<string, string> p1(v1.symbol, v1.formula), p2(v2.symbol, v2.formula);
    return p1 < p2;
}

struct InstanceAtom1stNeighboursNeighbours {
    int atomIdx, structureIdx;
    vector<Atoms1stNeighbours> neighbors;
};

bool less2(const InstanceAtom1stNeighboursNeighbours &v1,
           const InstanceAtom1stNeighboursNeighbours &v2) {
    return v1.neighbors < v2.neighbors;
}

void printDetailedStats(
    const string &fName, const string &header,
    const vector<pair<string, vector<int> > > &multitypes,
    const vector<vector<vector<pair<int, int> > > > &typeInstances,
    const vector<string> &structureNames,
    const vector<StructureWithDescriptors> structuresDescriptors,
    int nAssignedAtoms, int nNotAssignedAtoms, int nConsideredStructures,
    int nCompletelyRecognizedStructures) {
    ofstream out(fName);
    int nTypes, typeIdx, i, nInstances, nSubtypes, subtypeIdx, atomIdx,
        structureIdx, nNeighbours, neighbourIdx;

    if (!out.good())
        on_error::throwException(
            string("cannot print to file '") + fName + string("'"),
            __FILE__,
            __LINE__);

    string hash80(80, '#');
    out << hash80 << "\n"
        << header << hash80 << "\n\n"
        << "   content:\n"
        << "      (1) list of types which instances were found\n"
        << "      (2) list of types which instances were not found\n"
        << "      (3) fraction and number of atoms with assigned atom type\n"
        << "      (4) fraction and number of structures with all atoms with "
           "atom type assigned\n"
        << "      (5) number of instances per type\n"
        << "      (6) list of instances for each type \n"
        << "\n\n"
        << hash80 << "\n\n\n";

    //------------------------------------

    vector<int> nInstancesPerMultitype;
    nTypes = multitypes.size();
    nInstancesPerMultitype.assign(nTypes, 0);
    for (typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        nSubtypes = typeInstances[typeIdx].size();
        for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++)
            nInstancesPerMultitype[typeIdx] +=
                typeInstances[typeIdx][subtypeIdx].size();
    }

    //------------------------------------

    out << "(1) list of types which instances were found\n\n";
    i = 0;
    for (typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        if (nInstancesPerMultitype[typeIdx] > 0) {
            out << multitypes[typeIdx].first << " ";
            i++;
            if ((i + 1) % 10 == 0) out << "\n";
        }
    }
    out << "\n";

    out << "(2) list of types which instances were not found\n\n";

    i = 0;
    for (typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        if (nInstancesPerMultitype[typeIdx] == 0) {
            out << multitypes[typeIdx].first << " ";
            i++;
            if ((i + 1) % 10 == 0) out << "\n";
        }
    }
    out << "\n";

    out << "\n(3) fraction and number of atoms with assigned atom type\n";
    double nAllAtoms = double(nAssignedAtoms + nNotAssignedAtoms);
    out << "    atom type assigned to " << setprecision(2) << fixed
        << nAssignedAtoms / nAllAtoms * 100.0 << "% of atoms ("
        << nAssignedAtoms << ")\n";

    out << "\n(4) fraction and number of structures with all atoms with atom "
           "type assigned\n";
    out << "    " << nCompletelyRecognizedStructures << " of "
        << nConsideredStructures << " with completely recognized atom types "
        << " - " << setprecision(2) << fixed
        << double(nCompletelyRecognizedStructures) /
               double(nConsideredStructures) * 100.0
        << "%\n";

    out << "\n(5) number of instances per type\n\n";

    for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
        out << multitypes[typeIdx].first << " " << setw(8)
            << nInstancesPerMultitype[typeIdx] << "\n";

    out << "(6) list of instances for each type\n\n";

    vector<InstanceAtom1stNeighboursNeighbours> instanceAtomsNeighbours;
    vector<vector<int> > shells;
    vector<int> neighborNeighbours;

    // {first neighbour symbols, it neighbours formula}
    pair<string, vector<string> > firstNeighborsFormula;
    set<int> secondNeighbors;
    string formula;
    vector<int> atomicNumbers;

    for (typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        out << "\n" << "TYPE " << multitypes[typeIdx].first << "\n";
        nInstances = nInstancesPerMultitype[typeIdx];
        out << " number of instances " << nInstances << "\n";

        // formula_of_2nd_neighbours_and_idx.resize(nInstances);
        instanceAtomsNeighbours.clear();
        instanceAtomsNeighbours.resize(nInstances);
        nSubtypes = typeInstances[typeIdx].size();
        int instancesCounter = 0;
        for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++) {
            nInstances = typeInstances[typeIdx][subtypeIdx].size();

            for (i = 0; i < nInstances; i++) {
                structureIdx = typeInstances[typeIdx][subtypeIdx][i].first;
                atomIdx = typeInstances[typeIdx][subtypeIdx][i].second;
                nNeighbours = structuresDescriptors[structureIdx]
                                  .connectivity[atomIdx]
                                  .size();

                instanceAtomsNeighbours[instancesCounter].atomIdx = atomIdx;
                instanceAtomsNeighbours[instancesCounter].structureIdx =
                    structureIdx;
                instanceAtomsNeighbours[instancesCounter].neighbors.resize(
                    nNeighbours);
                firstNeighborsFormula.second.clear();
                for (int j = 0; j < nNeighbours;
                     j++)  // central atom 1-st neighbors
                {
                    neighbourIdx = structuresDescriptors[structureIdx]
                                       .connectivity[atomIdx][j];
                    instanceAtomsNeighbours[instancesCounter]
                        .neighbors[j]
                        .symbol = firstNeighborsFormula.first =
                        periodic_table::symbol(
                            structuresDescriptors[structureIdx]
                                .atomDescriptors[neighbourIdx]
                                .atomicNumber);
                    const vector<int> &neighbours =
                        structuresDescriptors[structureIdx]
                            .connectivity[neighbourIdx];
                    atomicNumbers.clear();
                    for (auto &neighbor :
                         neighbours)  // 1-st neighbor neighbors
                        atomicNumbers.push_back(
                            structuresDescriptors[structureIdx]
                                .atomDescriptors[neighbor]
                                .atomicNumber);

                    instanceAtomsNeighbours[instancesCounter]
                        .neighbors[j]
                        .formula = atomicNumbers2Formula(atomicNumbers);
                }
                sort(
                    instanceAtomsNeighbours[instancesCounter].neighbors.begin(),
                    instanceAtomsNeighbours[instancesCounter].neighbors.end(),
                    less1);
                instancesCounter++;
            }
        }

        std::sort(instanceAtomsNeighbours.begin(),
                  instanceAtomsNeighbours.end(),
                  less2);

        nInstances = instanceAtomsNeighbours.size();
        // formula, instance
        map<string, vector<pair<int, int> > > formulaGrouppedInstances;
        for (i = 0; i < nInstances; i++) {
            formula.clear();
            // ----- 1-st neighbors neighbors formula as string
            for (int neighIdx = 0;
                 neighIdx < instanceAtomsNeighbours[i].neighbors.size();
                 neighIdx++) {
                if (neighIdx != 0) formula += ";";
                formula +=
                    instanceAtomsNeighbours[i].neighbors[neighIdx].symbol +
                    string("(") +
                    instanceAtomsNeighbours[i].neighbors[neighIdx].formula +
                    string(")");
            }

            // -------
            structureIdx = instanceAtomsNeighbours[i].structureIdx;
            atomIdx = instanceAtomsNeighbours[i].atomIdx;
            formulaGrouppedInstances[formula].push_back(
                {structureIdx, atomIdx});
        }
        vector<pair<int, string> > neighbourGroups;
        for (auto &formula : formulaGrouppedInstances)
            neighbourGroups.push_back({formula.second.size(), formula.first});
        sort(neighbourGroups.begin(),
             neighbourGroups.end(),
             greater<pair<int, string> >());

        for (auto const &group : neighbourGroups) {
            const vector<pair<int, int> > &cases =
                formulaGrouppedInstances[group.second];
            out << "  " << group.first
                << " instances with neighbours formula: " << group.second
                << "\n";

            for (auto &c : cases) {
                string label =
                    structureNames[c.first] + string(",") +
                    structuresDescriptors[c.first]
                        .atomDescriptors[c.second]
                        .label;  // stats[typeIdx].atomLabels[index].first +
                                 // string(",") +
                                 // stats[typeIdx].atomLabels[index].second;
                out << setw(label.size() + 6) << label << "\n";
            }
        }
    }

    out.close();
}

void printMoreDetailedDetailedStats(
    const string &fName, const string &header,
    const vector<pair<string, vector<int> > > &multitypes,
    const vector<vector<vector<pair<int, int> > > > &typeInstances,
    const vector<string> &structureNames,
    const vector<StructureWithDescriptors> structuresDescriptors) {
    ofstream out(fName);
    int nTypes, typeIdx, i, nInstances, nSubtypes, subtypeIdx, atomIdx,
        structureIdx, nNeighbours, neighbourIdx;

    if (!out.good())
        on_error::throwException(
            string("cannot print to file '") + fName + string("'"),
            __FILE__,
            __LINE__);

    string hash80(80, '#');
    out << hash80 << "\n"
        << header << hash80 << "\n\n"
        << "   content: list of instances for each type\n\n";

    //------------------------------------

    vector<int> nInstancesPerMultitype;
    nTypes = multitypes.size();
    nInstancesPerMultitype.assign(nTypes, 0);
    for (typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        nSubtypes = typeInstances[typeIdx].size();
        for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++)
            nInstancesPerMultitype[typeIdx] +=
                typeInstances[typeIdx][subtypeIdx].size();
    }

    //------------------------------------

    vector<InstanceAtom1stNeighboursNeighbours> instanceAtomsNeighbours;
    vector<vector<int> > shells;
    vector<int> neighborNeighbours;

    // {first neighbour symbols, it neighbours formula}
    pair<string, vector<string> > firstNeighborsFormula;
    set<int> secondNeighbors;
    string formula;
    vector<int> atomicNumbers;

    for (typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        out << "\n" << "TYPE " << multitypes[typeIdx].first << "\n";
        nInstances = nInstancesPerMultitype[typeIdx];
        out << " number of instances " << nInstances << "\n";

        // formula_of_2nd_neighbours_and_idx.resize(nInstances);
        instanceAtomsNeighbours.clear();
        instanceAtomsNeighbours.resize(nInstances);
        nSubtypes = typeInstances[typeIdx].size();
        int instancesCounter = 0;
        for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++) {
            nInstances = typeInstances[typeIdx][subtypeIdx].size();

            for (i = 0; i < nInstances; i++) {
                structureIdx = typeInstances[typeIdx][subtypeIdx][i].first;
                atomIdx = typeInstances[typeIdx][subtypeIdx][i].second;
                nNeighbours = structuresDescriptors[structureIdx]
                                  .connectivity[atomIdx]
                                  .size();
                // instanceAtomsNeighbours[i].instanceIdx = i;
                instanceAtomsNeighbours[instancesCounter].atomIdx = atomIdx;
                instanceAtomsNeighbours[instancesCounter].structureIdx =
                    structureIdx;
                instanceAtomsNeighbours[instancesCounter].neighbors.resize(
                    nNeighbours);
                firstNeighborsFormula.second.clear();
                for (int j = 0; j < nNeighbours;
                     j++)  // central atom 1-st neighbors
                {
                    neighbourIdx = structuresDescriptors[structureIdx]
                                       .connectivity[atomIdx][j];
                    instanceAtomsNeighbours[instancesCounter]
                        .neighbors[j]
                        .symbol = firstNeighborsFormula.first =
                        periodic_table::symbol(
                            structuresDescriptors[structureIdx]
                                .atomDescriptors[neighbourIdx]
                                .atomicNumber);
                    if (structuresDescriptors[structureIdx]
                            .atomDescriptors[neighbourIdx]
                            .planar == Tribool::True)
                        instanceAtomsNeighbours[instancesCounter]
                            .neighbors[j]
                            .symbol += ",planar";
                    if (structuresDescriptors[structureIdx]
                            .atomDescriptors[neighbourIdx]
                            .planar == Tribool::False)
                        instanceAtomsNeighbours[instancesCounter]
                            .neighbors[j]
                            .symbol += ",not planar";

                    vector<int> const &ringSizes =
                        structuresDescriptors[structureIdx]
                            .atomDescriptors[neighbourIdx]
                            .planarRingsSizes;

                    if (!ringSizes.empty()) {
                        instanceAtomsNeighbours[instancesCounter]
                            .neighbors[j]
                            .symbol += ",rings(";
                        for (int k = 0; k < ringSizes.size(); k++) {
                            if (k > 0)
                                instanceAtomsNeighbours[instancesCounter]
                                    .neighbors[j]
                                    .symbol += ",";
                            instanceAtomsNeighbours[instancesCounter]
                                .neighbors[j]
                                .symbol += to_string(ringSizes[k]);
                        }
                        instanceAtomsNeighbours[instancesCounter]
                            .neighbors[j]
                            .symbol += ")";
                    }

                    int n3Rings = structuresDescriptors[structureIdx]
                                      .atomDescriptors[neighbourIdx]
                                      .n3memberRings;
                    int n4Rings = structuresDescriptors[structureIdx]
                                      .atomDescriptors[neighbourIdx]
                                      .n4memberRings;

                    if (n3Rings != 0)
                        instanceAtomsNeighbours[instancesCounter]
                            .neighbors[j]
                            .symbol +=
                            string(",3 rings: ") + to_string(n3Rings);
                    if (n4Rings != 0)
                        instanceAtomsNeighbours[instancesCounter]
                            .neighbors[j]
                            .symbol +=
                            string(",4 rings: ") + to_string(n4Rings);

                    const vector<int> &neighbours =
                        structuresDescriptors[structureIdx]
                            .connectivity[neighbourIdx];
                    atomicNumbers.clear();
                    for (auto &neighbor :
                         neighbours)  // 1-st neighbor neighbors
                        atomicNumbers.push_back(
                            structuresDescriptors[structureIdx]
                                .atomDescriptors[neighbor]
                                .atomicNumber);

                    instanceAtomsNeighbours[instancesCounter]
                        .neighbors[j]
                        .formula = atomicNumbers2Formula(atomicNumbers);
                }
                sort(
                    instanceAtomsNeighbours[instancesCounter].neighbors.begin(),
                    instanceAtomsNeighbours[instancesCounter].neighbors.end(),
                    less1);
                instancesCounter++;
            }
        }

        std::sort(instanceAtomsNeighbours.begin(),
                  instanceAtomsNeighbours.end(),
                  less2);

        nInstances = instanceAtomsNeighbours.size();
        // formula, instance
        map<string, vector<pair<int, int> > > formulaGrouppedInstances;
        for (i = 0; i < nInstances; i++) {
            formula.clear();
            // ----- 1-st neighbors neighbors formula as string
            for (int neighIdx = 0;
                 neighIdx < instanceAtomsNeighbours[i].neighbors.size();
                 neighIdx++) {
                if (neighIdx != 0) formula += ";";
                formula +=
                    instanceAtomsNeighbours[i].neighbors[neighIdx].symbol +
                    string(", (") +
                    instanceAtomsNeighbours[i].neighbors[neighIdx].formula +
                    string(")");
            }

            // -------
            structureIdx = instanceAtomsNeighbours[i].structureIdx;
            atomIdx = instanceAtomsNeighbours[i].atomIdx;
            formulaGrouppedInstances[formula].push_back(
                {structureIdx, atomIdx});
        }
        vector<pair<int, string> > neighbourGroups;
        for (auto &formula : formulaGrouppedInstances)
            neighbourGroups.push_back({formula.second.size(), formula.first});
        sort(neighbourGroups.begin(),
             neighbourGroups.end(),
             greater<pair<int, string> >());

        for (auto const &group : neighbourGroups) {
            const vector<pair<int, int> > &cases =
                formulaGrouppedInstances[group.second];
            out << "  " << group.first
                << " instances with neighbours formula: " << group.second
                << "\n";

            for (auto &c : cases) {
                string label =
                    structureNames[c.first] + string(",") +
                    structuresDescriptors[c.first]
                        .atomDescriptors[c.second]
                        .label;  // stats[typeIdx].atomLabels[index].first +
                                 // string(",") +
                                 // stats[typeIdx].atomLabels[index].second;
                out << setw(label.size() + 6) << label << "\n";
            }
        }
    }

    out.close();
}

string secondNeighboursLabel(
    const StructureWithDescriptors &structureDescriptors, int &atomIdx) {
    vector<vector<int> > neighbours;
    graph_algorithms::breadth_first_search(
        structureDescriptors.connectivity, atomIdx, neighbours, 2);
    if (neighbours.size() < 3) return string("0");
    int n = neighbours[2].size();
    set<string> symbols;
    for (int i = 0; i < n; i++)
        symbols.insert(periodic_table::symbol(
            structureDescriptors.atomDescriptors[neighbours[2][i]]
                .atomicNumber));
    string result;
    for (auto &symbol : symbols) result += symbol;
    result += to_string(n);

    return result;
}

typedef tuple<int, string, string, string, int, int, string>
    UnassignedAtomDescriptors;

void printAssignmentInfo(
    const string &fName, const vector<string> &structureIds,
    const std::string &header, const DescriptorsSettings &descriptorsSettings,
    const vector<StructureWithDescriptors> &structureDescriptors,
    const vector<vector<int> > &typeIndices, const vector<AtomType> &types,
    // const vector<vector<LocalCoordinateSystem<int> > > &lcs,
    const vector<vector<string> > &lcs,
    map<UnassignedAtomDescriptors, vector<string> > &sortedUnassignedAtoms,
    const vector<int> &nAtomsInAsymmetricUnit) {
    sortedUnassignedAtoms.clear();
    UnassignedAtomDescriptors unassignedAtomDescriptors;
    ofstream out(fName);
    char xyz[] = {'X', 'Y', 'Z'};
    if (!out.good())
        on_error::throwException(
            string("cannot write to log file '") + fName + string("'"),
            __FILE__,
            __LINE__);

    int nStructures, nAtoms, atomIdx, structureIdx;
    string ringsSizesStr;
    vector<int> ringsSizes;
    string hash80(80, '#');
    nStructures = structureDescriptors.size();

    out << hash80 << "\n"
        << header << hash80 << "\n\n"
        << "       Structural Descriptors for each structure and atom and\n"
        << "       atom type assignemnt and local cordinate system "
           "information\n\n"
        << hash80 << "\n\n"
        << "       Number of structures " << nStructures << "\n\n"
        << hash80 << "\n\n";

    vector<string> labels;

    for (structureIdx = 0; structureIdx < nStructures; structureIdx++) {
        labels.clear();
        for (auto const &atomData :
             structureDescriptors[structureIdx].atomDescriptors)
            labels.push_back(atomData.label);

        out << "\n Structure " << structureIdx + 1 << " "
            << structureIds[structureIdx] << endl;
        out << "\n number of planar rings: "
            << structureDescriptors[structureIdx].planarRings.size() << endl;
        for (int ringIdx = 0;
             ringIdx < structureDescriptors[structureIdx].planarRings.size();
             ringIdx++) {
            out << "\n  ring " << ringIdx + 1 << "\n    planarity esd "
                << structureDescriptors[structureIdx]
                       .planarRingsPlanarityEsd[ringIdx];
            out << "\n    atoms: ";
            for (int i = 0;
                 i <
                 structureDescriptors[structureIdx].planarRings[ringIdx].size();
                 i++)
                out << " "
                    << labels[structureDescriptors[structureIdx]
                                  .planarRings[ringIdx][i]];
            out << "\n";
        }

        //------------------------------------------------------------------------------------

        nAtoms = nAtomsInAsymmetricUnit[structureIdx];
        out << "\n Atomic Descriptors\n";
        out << "\n number of atoms " << nAtoms << endl;

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
            string symbol;
            vector<int> _v{structureDescriptors[structureIdx]
                               .atomDescriptors[atomIdx]
                               .neighborsFormula.begin(),
                           structureDescriptors[structureIdx]
                               .atomDescriptors[atomIdx]
                               .neighborsFormula.end()};
            string formulaStr = atomicNumbers2Formula(_v);
            string planarityStr;

            planarityStr = atom_typing_utilities::planarity_as_string(
                structureDescriptors[structureIdx]
                    .atomDescriptors[atomIdx]
                    .planar);

            out << setw(8) << labels[atomIdx] << " atomic number " << setw(4)
                << structureDescriptors[structureIdx]
                       .atomDescriptors[atomIdx]
                       .atomicNumber
                << " , neighbours formula " << setw(12) << formulaStr
                << " , is planar " << setw(15) << planarityStr
                << " , planarity esd " << setw(10) << setprecision(6) << fixed
                << structureDescriptors[structureIdx]
                       .atomDescriptors[atomIdx]
                       .planarityDisplacementEsd
                << " , belongs to planar rings: ";

            ringsSizesStr.clear();

            if (structureDescriptors[structureIdx]
                    .atomDescriptors[atomIdx]
                    .planarRingsIndices.empty()) {
                out << "  no ring size: ";
                ringsSizesStr = " - ";
            } else {
                out << " yes ring size: ";
                ringsSizes.clear();
                const vector<int> &ringsIndices =
                    structureDescriptors[structureIdx]
                        .atomDescriptors[atomIdx]
                        .planarRingsIndices;
                for (auto ringIdx : structureDescriptors[structureIdx]
                                        .atomDescriptors[atomIdx]
                                        .planarRingsIndices)
                    ringsSizes.push_back(structureDescriptors[structureIdx]
                                             .planarRings[ringIdx]
                                             .size());
                sort(ringsSizes.begin(), ringsSizes.end());
                for (auto ringSize : ringsSizes) {
                    if (!ringsSizesStr.empty()) ringsSizesStr += ",";
                    ringsSizesStr += to_string(ringSize);
                }
            }

            out << setw(11) << ringsSizesStr;

            // 3 and 4 member rings

            out << " , in "
                << structureDescriptors[structureIdx]
                       .atomDescriptors[atomIdx]
                       .n3memberRings
                << " 3-member rings "
                << " , in "
                << structureDescriptors[structureIdx]
                       .atomDescriptors[atomIdx]
                       .n4memberRings
                << " 4-member rings ";

            // collect information on unassigned atoms

            if (typeIndices[structureIdx][atomIdx] < 0) {
                // atomic number, neighbours, planarity, rings
                std::get<0>(unassignedAtomDescriptors) =
                    structureDescriptors[structureIdx]
                        .atomDescriptors[atomIdx]
                        .atomicNumber;
                std::get<1>(unassignedAtomDescriptors) = formulaStr;
                std::get<2>(unassignedAtomDescriptors) = planarityStr;
                std::get<3>(unassignedAtomDescriptors) = ringsSizesStr;
                std::get<4>(unassignedAtomDescriptors) =
                    structureDescriptors[structureIdx]
                        .atomDescriptors[atomIdx]
                        .n3memberRings;
                std::get<5>(unassignedAtomDescriptors) =
                    structureDescriptors[structureIdx]
                        .atomDescriptors[atomIdx]
                        .n4memberRings;
                std::get<6>(unassignedAtomDescriptors) = secondNeighboursLabel(
                    structureDescriptors[structureIdx], atomIdx);

                string atomLabel =
                    structureIds[structureIdx] + string(",") + labels[atomIdx];
                sortedUnassignedAtoms[unassignedAtomDescriptors].push_back(
                    atomLabel);
            }

            // ---------------------------------------

            out << " , neighbours ";
            const vector<int> &neighbourIndices =
                structureDescriptors[structureIdx].connectivity[atomIdx];
            vector<string> neighborsLabels;

            for (int nIdx = 0; nIdx < neighbourIndices.size(); nIdx++)
                neighborsLabels.push_back(
                    structureDescriptors[structureIdx]
                        .atomDescriptors[neighbourIndices[nIdx]]
                        .label);
            sort(neighborsLabels.begin(), neighborsLabels.end());

            for (int nIdx = 0; nIdx < neighborsLabels.size(); nIdx++) {
                if (nIdx > 0) out << ",";
                out << neighborsLabels[nIdx];
            }

            vector<string> neighborsTypes;
            out << " , neighbours types ";

            for (int nIdx = 0; nIdx < neighbourIndices.size(); nIdx++) {
                int neighbourIdx = neighbourIndices[nIdx];
                string typeName;
                if (neighbourIdx >= nAtomsInAsymmetricUnit[structureIdx])
                    typeName = "?";
                else
                    typeName =
                        (typeIndices[structureIdx][neighbourIdx] >= 0
                             ? types[typeIndices[structureIdx][neighbourIdx]].id
                             : string("-"));
                neighborsTypes.push_back(typeName);
            }

            sort(neighborsTypes.begin(), neighborsTypes.end());

            for (int nIdx = 0; nIdx < neighbourIndices.size(); nIdx++) {
                if (nIdx > 0) out << ",";
                out << neighborsTypes[nIdx];
            }

            out << endl;
        }  // for(int atomIdx = 0; atomIdx <
           // structreDescriptors[moleculeIdx].atomDescriptors.size();
           // atomIdx++)

        out << "\n\n Atom Type Assignemnt\n\n";
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
            string typeName =
                (typeIndices[structureIdx][atomIdx] >= 0
                     ? types[typeIndices[structureIdx][atomIdx]].id
                     : string("-"));
            out << setw(8) << " " << labels[atomIdx] << " " << setw(8)
                << typeName << " ";

            // print lcs
            if (typeIndices[structureIdx][atomIdx] >= 0)
                out << lcs
                        [structureIdx]
                        [atomIdx];  // ubdbLcsAsString(lcs[structureIdx][atomIdx],
                                    // labels);
            out << endl;
        }

    }  // for (moleculeIdx = 0; moleculeIdx < nMolecules; moleculeIdx++)
}

void printSettingsInfo(const string &fName, const string &header,
                       const BankSettings &settings) {
    ofstream out(fName);

    if (!out.good())
        on_error::throwException(
            string("cannot print to file '") + fName + string("'"),
            __FILE__,
            __LINE__);

    string hash80(80, '#');
    hash80 += string("\n");

    out << hash80 << header << hash80 << "\n";

    out << "the following settings were used:"
        << "\n      ring considered not planar if it contains atoms with more "
           "than "
        << settings.descriptorsSettings.atomInRingMaxNeighbourCount
        << " neighbours"
        << "\n      ring considered not planar if it contains atoms with "
           "planarity above "
        << settings.descriptorsSettings.atomInRingPlanarityThreshold
        << "\n      ring considered not planar if its planarity is above "
        << settings.descriptorsSettings.ringPlanarityThreshold
        << "\n      atom considered planar if the planarity is below "
        << settings.descriptorsSettings.atomPlanarityThreshold
        << "\n      atoms are bonded if the interatomic distance is below sum "
           "of their covalent bonds plus "
        << settings.descriptorsSettings.covalentBondThreshold
        << "\n      minimal absolute value of P_lm to be included in the bank "
           "is "
        << settings.min_plm
        << "\n      P_lm is included in the bank if its absolute value is at "
           "least "
        << settings.nSigma << " times higher than the P_lm standard deviation"
        << "\n      minimal number of atom type instances for parameters "
           "update "
        << settings.min_n_instaces << "\n\n"
        << hash80 << "\n";
    out.close();
}

void printUnassignedAtomsInfo(
    const string fName, string header,
    map<UnassignedAtomDescriptors, vector<string> > &sortedUnassignedAtoms) {
    ofstream out(fName);
    if (!out.good())
        on_error::throwException(
            string("cannot write to log file '") + fName + string("'"),
            __FILE__,
            __LINE__);

    string hash80(80, '#');

    out << hash80 << "\n"
        << header << hash80 << "\n\n"
        << "       atoms with unassigned types sorted by descriptors\n\n"
        << hash80 << "\n\n";
    int i, n;
    for (auto &atoms : sortedUnassignedAtoms) {
        out << "atomic number " << setw(4) << get<0>(atoms.first)
            << " neighbours " << setw(12) << get<1>(atoms.first)
            << " planarity " << setw(15) << get<2>(atoms.first)
            << " planar rings with planar atoms " << setw(15)
            << get<3>(atoms.first) << " ,3 member rings " << setw(3)
            << get<4>(atoms.first) << " ,4 member rings " << setw(3)
            << get<5>(atoms.first) << " , 2-nd neighbours "
            << get<6>(atoms.first) << "\n\n";

        n = atoms.second.size();
        out << " number of atoms : " << n << "\n";
        for (i = 0; i < n; i++) {
            out << atoms.second[i] << " ";
            if ((i + 1) % 8 == 0) out << "\n";
        }
        if (i % 8 != 0) out << "\n";
        out << "\n-------------------------------------------------------------"
               "--------------\n\n";
    }

    out.close();
}

void printUnassignedAtomsInfoSortedByN(
    const string fName, string header,
    map<UnassignedAtomDescriptors, vector<string> > &sortedUnassignedAtoms) {
    ofstream out(fName);
    if (!out.good())
        on_error::throwException(
            string("cannot write to log file '") + fName + string("'"),
            __FILE__,
            __LINE__);

    string hash80(80, '#');

    out << hash80 << "\n"
        << header << hash80 << "\n\n"
        << "       atoms with unassigned types sorted by number of atoms in "
           "group\n\n"
        << hash80 << "\n\n";

    vector<pair<UnassignedAtomDescriptors, vector<string> > >
        unassignedAtomGroups;
    vector<pair<int, int> > groupMultiplicityAndIndex;

    int index = 0;
    for (auto &atoms : sortedUnassignedAtoms) {
        unassignedAtomGroups.push_back(atoms);
        groupMultiplicityAndIndex.push_back({atoms.second.size(), index});
        index++;
    }

    sort(groupMultiplicityAndIndex.begin(),
         groupMultiplicityAndIndex.end(),
         greater<pair<int, int> >());

    for (int i = 0; i < groupMultiplicityAndIndex.size(); i++) {
        index = groupMultiplicityAndIndex[i].second;

        out << "atomic number " << setw(4)
            << get<0>(unassignedAtomGroups[index].first) << " neighbours "
            << setw(12) << get<1>(unassignedAtomGroups[index].first)
            << " planarity " << setw(15)
            << get<2>(unassignedAtomGroups[index].first)
            << " planar rings with planar atoms " << setw(15)
            << get<3>(unassignedAtomGroups[index].first) << " ,3 member rings "
            << setw(3) << get<4>(unassignedAtomGroups[index].first)
            << " ,4 member rings " << setw(3)
            << get<5>(unassignedAtomGroups[index].first)
            << " , 2-nd neighbours "
            << get<6>(unassignedAtomGroups[index].first) << "\n\n";

        int n = unassignedAtomGroups[index].second.size();
        out << " number of atoms : " << n << "\n";
        int j;
        for (j = 0; j < n; j++) {
            out << setw(20) << unassignedAtomGroups[index].second[j] << " ";
            if ((j + 1) % 8 == 0) out << "\n";
        }
        if (j % 8 != 0) out << "\n";
        out << "\n-------------------------------------------------------------"
               "--------------\n\n";
    }

    out.close();
}

bool read_bank_and_assign_atoms(
    const string bank_filepath, const Crystal &crystal,
    vector<AtomType> &atomTypes, vector<AtomTypeHC_Parameters> &hcParameters,
    BankSettings &bankSettings, CrystalAtomTypeAssigner &crystalAssigner,
    int &nAtomsInAsymmetricUnit, StructureWithDescriptors &structure,
    vector<int> &typeIds, vector<LocalCoordinateSystem<AtomInCrystalID> > &lcs,
    vector<string> &lcs_strings) {
    // read bank parameters
    MATTS_BankReader bankReader;

    ifstream bankStream{bank_filepath};
    bankReader.read(bankStream, atomTypes, hcParameters, bankSettings, true);

    crystalAssigner.setDescriptorsSettings(bankSettings.descriptorsSettings);
    crystalAssigner.setAtomTypes(atomTypes);

    vector<int> atomicNumbers;
    vector<string> labels;
    bool processingCrystalStructure = true;

    nAtomsInAsymmetricUnit = crystal.atoms.size();

    crystalAssigner.assign(crystal, typeIds, lcs, structure);

    for (auto &atom : crystal.atoms) labels.push_back(atom.label);

    for (auto &coordinateSystem : lcs)
        lcs_strings.push_back(ubdbLcsAsString(coordinateSystem, labels));

    return true;
}

void findMultitypes(const vector<AtomType> &types,
                    vector<pair<string, vector<int> > > &multitypes) {
    multitypes.clear();
    map<string, vector<int> > mtypes;

    for (int i = 0; i < types.size(); i++) mtypes[types[i].id].push_back(i);

    for (auto mtype : mtypes) multitypes.push_back({mtype.first, mtype.second});
}

void findTypeInstances(const vector<pair<string, vector<int> > > &multitypes,
                       const vector<vector<int> > &typeIndices,
                       vector<vector<vector<pair<int, int> > > > &typeInstances,
                       int total_n_subtypes) {
    int multitypeIdx, subtypeIdx, nSubtypes, nMultiTypes = multitypes.size();
    typeInstances.clear();
    typeInstances.resize(nMultiTypes);
    for (multitypeIdx = 0; multitypeIdx < nMultiTypes; multitypeIdx++)
        typeInstances[multitypeIdx].resize(
            multitypes[multitypeIdx].second.size());

    vector<pair<int, int> > type2multitypeSubtype(total_n_subtypes);

    for (multitypeIdx = 0; multitypeIdx < nMultiTypes; multitypeIdx++) {
        nSubtypes = multitypes[multitypeIdx].second.size();

        for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++)
            type2multitypeSubtype[multitypes[multitypeIdx].second[subtypeIdx]] =
                {multitypeIdx, subtypeIdx};
    }
    //
    int structureIdx, nStructures, atomIdx, nAtoms;

    nStructures = typeIndices.size();

    for (structureIdx = 0; structureIdx < nStructures; structureIdx++) {
        nAtoms = typeIndices[structureIdx].size();
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if (typeIndices[structureIdx][atomIdx] >= 0) {
                multitypeIdx =
                    type2multitypeSubtype[typeIndices[structureIdx][atomIdx]]
                        .first;
                subtypeIdx =
                    type2multitypeSubtype[typeIndices[structureIdx][atomIdx]]
                        .second;
                typeInstances[multitypeIdx][subtypeIdx].push_back(
                    {structureIdx, atomIdx});
            }
    }
}

void write_assignment_logs(
    const Crystal crystal, const vector<AtomType> atomTypes,
    const vector<AtomTypeHC_Parameters> hcParameters,
    const BankSettings bankSettings,
    const CrystalAtomTypeAssigner crystalAssigner,
    const int nAtomsInAsymmetricUnit, const StructureWithDescriptors structure,
    const vector<int> typeIds,
    const vector<LocalCoordinateSystem<AtomInCrystalID> > lcs,
    const vector<string> lcs_strings) {
    int nConsideredStructures, nCompletelyRecognizedStructures;

    try {
        ofstream default_log{"default_assignment_log.txt"};
        crystalAssigner.printAssignment(default_log, crystal, typeIds, lcs);
        default_log.close();

        vector<string> structureIds{"Structure"};
        vector<vector<int> > typeIndices;
        vector<int> natomsInAsymmetricUnit;
        vector<vector<string> > lcs;

        int structureIdx, nStructures = 1;

        vector<StructureWithDescriptors> structureDescriptors;

        vector<string> validStructureIds;
        nConsideredStructures = nCompletelyRecognizedStructures = 0;

        natomsInAsymmetricUnit.push_back(nAtomsInAsymmetricUnit);
        validStructureIds.push_back(structureIds[0]);
        structureDescriptors.push_back(structure);
        lcs.push_back(lcs_strings);
        typeIndices.push_back(typeIds);
        nConsideredStructures++;

        nCompletelyRecognizedStructures =
            (find(typeIds.begin(), typeIds.end(), -1) == typeIds.end());

        validStructureIds.swap(structureIds);

        int nAssignedAtoms, nNotAssignedAtoms;
        nAssignedAtoms = nNotAssignedAtoms = 0;
        for (auto const &atomTypeIndices : typeIndices)
            for (int typeIdx : atomTypeIndices)
                if (typeIdx < 0)
                    nNotAssignedAtoms++;
                else
                    nAssignedAtoms++;

        vector<vector<vector<pair<int, int> > > > typeInstances;
        vector<pair<string, vector<int> > >
            multitypes;  //[i].first - label [i].second - list of indices of
                         // types

        findMultitypes(atomTypes, multitypes);
        findTypeInstances(
            multitypes, typeIndices, typeInstances, atomTypes.size());

        printDetailedStats("type_instances_info.txt",
                           "",
                           multitypes,
                           typeInstances,
                           structureIds,
                           structureDescriptors,
                           nAssignedAtoms,
                           nNotAssignedAtoms,
                           nConsideredStructures,
                           nCompletelyRecognizedStructures);

        printMoreDetailedDetailedStats("type_instances_detailed_info.txt",
                                       "",
                                       multitypes,
                                       typeInstances,
                                       structureIds,
                                       structureDescriptors);

        map<UnassignedAtomDescriptors, vector<string> > sortedUnassignedAtoms;

        printAssignmentInfo("assignment_info.txt",
                            structureIds,
                            "",
                            bankSettings.descriptorsSettings,
                            structureDescriptors,
                            typeIndices,
                            atomTypes,
                            lcs,
                            sortedUnassignedAtoms,
                            natomsInAsymmetricUnit);

        printUnassignedAtomsInfo(
            "unassigned_atoms_info.txt", "", sortedUnassignedAtoms);
        printUnassignedAtomsInfoSortedByN(
            "unassigned_atoms_sorted_by_n.txt", "", sortedUnassignedAtoms);
        printSettingsInfo("settings_info.txt", "", bankSettings);
    } catch (exception &e) {
        cout << e.what() << endl;
    }
}

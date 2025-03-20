#pragma once

#include <string>
#include <vector>

#include "discamb/AtomTyping/AtomType.h"
#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystem.h"
#include "discamb/AtomTyping/StructureWithDescriptors.h"
#include "discamb/CrystalStructure/AtomInCrystalID.h"
#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/IO/MATTS_BankReader.h"

bool read_bank_and_assign_atoms(
    const std::string bank_filepath, const discamb::Crystal &crystal,
    std::vector<discamb::AtomType> &atomTypes,
    std::vector<discamb::AtomTypeHC_Parameters> &hcParameters,
    discamb::BankSettings &bankSettings,
    discamb::CrystalAtomTypeAssigner &crystalAssigner,
    int &nAtomsInAsymmetricUnit, discamb::StructureWithDescriptors &structure,
    std::vector<int> &typeIds,
    std::vector<discamb::LocalCoordinateSystem<discamb::AtomInCrystalID> > &lcs,
    std::vector<std::string> &lcs_strings);

void write_assignment_logs(
    const discamb::Crystal crystal,
    const std::vector<discamb::AtomType> atomTypes,
    const std::vector<discamb::AtomTypeHC_Parameters> hcParameters,
    const discamb::BankSettings bankSettings,
    const discamb::CrystalAtomTypeAssigner crystalAssigner,
    const int nAtomsInAsymmetricUnit,
    const discamb::StructureWithDescriptors structure,
    const std::vector<int> typeIds,
    const std::vector<discamb::LocalCoordinateSystem<discamb::AtomInCrystalID> >
        lcs,
    const std::vector<std::string> lcs_strings);

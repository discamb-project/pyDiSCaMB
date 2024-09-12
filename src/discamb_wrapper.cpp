#include "discamb_wrapper.hpp"

#include "discamb/BasicUtilities/discamb_version.h"

std::string get_discamb_version(){

    return "DiSCaMB version: " + discamb::discamb_version::version();
} 

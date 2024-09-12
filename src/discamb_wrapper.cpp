#include "discamb_wrapper.hpp"

#include <iostream>

int main(){
    std::cout << "DiSCaMB version: " << discamb::discamb_version::version() << std::endl;
    return 0;
}

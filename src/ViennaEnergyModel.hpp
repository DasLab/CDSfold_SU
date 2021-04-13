#include "EnergyModel.hpp"

#pragma once

#ifdef USE_VIENNA_ENERGY_MODEL

class ViennaEnergyModel: public EnergyModel {

public:

    /* constructor - initialize the energy parameters here */
    ViennaEnergyModel() : EnergyModel((scaledEnergyParams*) scale_parameters()) {};
   
    void repr() override {
        std::cout << "Vienna Energy Model" << std::endl;
    }
  
    void updateEnergyFoldParams() override {
        update_fold_params();
    }
    
    int TermAU(int const &type);
    int E_hairpin(int size, int type, int si1, int sj1, const char *string);
    int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1);
};

#endif

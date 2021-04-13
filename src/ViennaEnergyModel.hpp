#include "EnergyModel.hpp"

#pragma once

extern "C" {
    #include "params.h"
    #include "utils.h"
    #include "energy_const.h"
    #include "fold.h"
}

class ViennaEnergyModel: public EnergyModel {

protected:
    vrna_param_s* energyParams_;

public:
    /* constructor - initialize the energy parameters here */
    ViennaEnergyModel() : energyParams_(scale_parameters()) {};
    
    void repr() override {
        std::cout << "Vienna Energy Model" << std::endl;
    }
  
    void updateEnergyFoldParams() override {
        update_fold_params();
    }
   
    /* energy functions */
    int TermAU(int const &type);
    int E_hairpin(int size, int type, int si1, int sj1, const char *string);
    int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1);

    /* accessor functions */
    int getTerminalAU() override { return energyParams_->TerminalAU;}
    int getMLbase() override { return energyParams_->MLbase;}
    int getMLintern(unsigned int i) override { return energyParams_->MLintern[i];}
    int getMLclosing() override { return energyParams_->MLclosing;}

};


#include "EnergyModel.hpp"

#pragma once

/*
 * struct for holding RNA energy parameters for the dummy model
 */

/* template dummy code to suppress argument not used warnings */
template <typename... Targs>
void DUMMY_CODE(Targs &&... /* unused */){}


class dummyEnergyParams : public energyParams {
public:
    int       MLbase;
    int       MLintern[NBPAIRS + 1];
    int       MLclosing;
    int       TerminalAU;
};


/* DummyEnergyModel
 *
 * Empty energy model that is another child of EnergyModel; sibling of ViennaEnergyModel 
 * used to demonstrate that CDSFold can be used without Vienna */
class DummyEnergyModel: public EnergyModel {

protected:
   dummyEnergyParams* energyParams_; 

public:
    DummyEnergyModel() {};

    void repr() override {
        std::cout << "Dummy Energy Model" << std::endl;
    }

    void updateEnergyFoldParams() override {}
    
    //dummyEnergyParams* getEnergyParams() override {
    //    return energyParams_;
    //}

    /* energy functions */
    int TermAU(int const &type) override {
        DUMMY_CODE(type); 
        return 0;
    }
    
    int E_hairpin(int size, int type, int si1, int sj1, const char *string) override {
        DUMMY_CODE(size, type, si1, sj1, string); 
        return 0;
    }
    
    int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1) override {
        DUMMY_CODE(n1, n2, type, type_2, si1, sj1, sp1, sq1); 
        return 0;
    }

    /* accessor functions */
    int getTerminalAU() override { return energyParams_->TerminalAU;}
    int getMLbase() override { return energyParams_->MLbase;}
    int getMLintern(unsigned int i) override { return energyParams_->MLintern[i];}
    int getMLclosing() override { return energyParams_->MLclosing;}

};


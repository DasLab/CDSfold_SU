/*
 * DummyEnergyModel.hpp
 *
 * Header file for a dummy energy model. This class is an alternative to 
 * the ViennaEnergyModel and demonstrates building and compiling cdsfold
 * without Vienna RNA. This energy model will not correctly backtrace
 */

#include "EnergyModel.hpp"

#pragma once


/* template dummy code to suppress argument not used warnings */
template <typename... Targs>
void DUMMY_CODE(Targs &&... /* unused */){}


/* minimum required energy parameters*/
class dummyEnergyParams {
public:
    int       MLbase;
    int       MLintern[NBPAIRS + 1];
    int       MLclosing;
    int       TerminalAU;
};


/* DummyEnergyModel
 *
 */

class DummyEnergyModel: public EnergyModel {

protected:
    /* pointer to energy parameters */
    std::unique_ptr<dummyEnergyParams> energyParams_; 

public:
    DummyEnergyModel() {};

    void repr() override {
        std::cout << "Dummy Energy Model" << std::endl;
    }

    void updateEnergyFoldParams() override {}

    /* dummy energy functions - all return 0 */
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


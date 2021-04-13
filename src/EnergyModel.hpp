/*
 * Class for the energy model
 *
 */
#pragma once

#include "constants.hpp"
#include <iostream>
#include <ostream>
#include <memory>


/* energyParams
 *
 * base class for energy parameters for other models. This class is empty
 * but is inherited by the energy parameters for other energy models. The
 * energyParams dummy class is needed so that the return type of getEnergyParams
 * in energyModel can be covariant. */

class energyParams {

};


/* base class representing an energy model */
class EnergyModel {

public:
    /* constructor */ 
    EnergyModel() {};
    
    /* repr function - prints name of the class */
    virtual void repr() = 0;

    /* updateEnergyFoldParams - updates the fold parameters */
    virtual void updateEnergyFoldParams() = 0;

    /* energy functions */
    virtual int TermAU(int const &type) = 0;
    virtual int E_hairpin(int size, int type, int si1, int sj1, const char *string) = 0;
    virtual int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1) = 0;
    
    /* accessor functions for certain parameters */
    virtual int getTerminalAU() = 0;
    virtual int getMLbase() = 0;
    virtual int getMLintern(unsigned int i) = 0;
    virtual int getMLclosing() = 0;

};



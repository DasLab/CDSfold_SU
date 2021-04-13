/*
 * Class for the energy model
 *
 */

#include "constants.hpp"
#include <iostream>
#include <ostream>
#include <memory>

/* Preprocessor directives for compiling an energy model
 */ 
#pragma once
#define USE_VIENNA_ENERGY_MODEL

#ifdef USE_VIENNA_ENERGY_MODEL

extern "C" {
    #include "params.h"
    #include "utils.h"
    #include "energy_const.h"
    #include "fold.h"
}

#endif


class EnergyParams {
public:
    int    TerminalAU;
    int    MLbase;
    int    MLintern[NBPAIRS + 1];
    int    MLclosing;
};


/* base class representing an energy model */
class EnergyModel {

public:
    /* constructor 
     *             into the energyParams_ attribute */ 
    EnergyModel() {};
    
    /* getEnergyParams - returns a shared pointer to the energy parameters
     */
    virtual EnergyParams* getEnergyParams() = 0; 
    /* repr function - prints name of the class */
    virtual void repr() = 0;

    /* updateEnergyFoldParams - updates the fold parameters */
    virtual void updateEnergyFoldParams() = 0;

    virtual int TermAU(int const &type) = 0;
    virtual int E_hairpin(int size, int type, int si1, int sj1, const char *string) = 0;
    virtual int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1) = 0;

};


/* Need something like vrna_param_s type struct for holding temperature scaled energy
 * parameters - typedef'd as paramT and attribute of Problem */

class DummyEnergyModel: public EnergyModel {
public:
    DummyEnergyModel();

    void repr() override {
        std::cout << "Dummy Energy Model" << std::endl;
    }

    void updateEnergyFoldParams() override {

    }

};


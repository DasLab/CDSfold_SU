#include "EnergyModel.hpp"

#pragma once

#ifdef USE_VIENNA_ENERGY_MODEL

/*
 * struct for holding RNA scaled energy parameters
 * from ViennaRNA/include/ViennaRNA/params/basic.h vrna_param_s
 */


class scaledEnergyParams : public energyParams{
public:
    int       id;
    int       stack[NBPAIRS + 1][NBPAIRS + 1];
    int       hairpin[31];
    int       bulge[MAXLOOP + 1];
    int       internal_loop[MAXLOOP + 1];
    int       mismatchExt[NBPAIRS + 1][5][5];
    int       mismatchI[NBPAIRS + 1][5][5];
    int       mismatch1nI[NBPAIRS + 1][5][5];
    int       mismatch23I[NBPAIRS + 1][5][5];
    int       mismatchH[NBPAIRS + 1][5][5];
    int       mismatchM[NBPAIRS + 1][5][5];
    int       dangle5[NBPAIRS + 1][5];
    int       dangle3[NBPAIRS + 1][5];
    int       int11[NBPAIRS + 1][NBPAIRS + 1][5][5];
    int       int21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
    int       int22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
    int       ninio[5];
    double    lxc;
    int       MLbase;
    int       MLintern[NBPAIRS + 1];
    int       MLclosing;
    int       TerminalAU;
    int       DuplexInit;
    int       Tetraloop_E[200];
    char      Tetraloops[1401];
    int       Triloop_E[40];
    char      Triloops[241];
    int       Hexaloop_E[40];
    char      Hexaloops[1801];
    int       TripleC;
    int       MultipleCA;
    int       MultipleCB;
    int       gquad [VRNA_GQUAD_MAX_STACK_SIZE + 1]
    [3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];

    double    temperature;      /**<  @brief  Temperature used for loop contribution scaling */

    vrna_md_t model_details;    /**<  @brief  Model details to be used in the recursions */
    char      param_file[256];  /**<  @brief  The filename the parameters were derived from, or empty string if they represent the default */

};

class ViennaEnergyModel: public EnergyModel {

protected:
    scaledEnergyParams* energyParams_;

public:

    /* constructor - initialize the energy parameters here */
    ViennaEnergyModel() : energyParams_((scaledEnergyParams*) scale_parameters()) {};
  
    scaledEnergyParams* getEnergyParams() override {
        return energyParams_;
    }
    
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

#endif

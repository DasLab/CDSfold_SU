/*
 * Class for the energy model
 *
 */


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


/*
 * struct for holding RNA scaled energy parameters
 * from ViennaRNA/include/ViennaRNA/params/basic.h vrna_param_s
 */

#define NBPAIRS 7
#define MAXLOOP 30
#define VRNA_GQUAD_MAX_STACK_SIZE 7
#define VRNA_GQUAD_MAX_LINKER_LENGTH 15


typedef struct scaledEnergyParams {
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

} scaledEnergyParams;


/* base class representing an energy model */
class EnergyModel {

public:
    unique_ptr<scaledEnergyParams> P_;  // TODO - make this protected
    EnergyModel(scaledEnergyParams* p = nullptr) : P_(p) {}
    //virtual void update_fold_params() = 0;

};


/* Need something like vrna_param_s type struct for holding temperature scaled energy
 * parameters - typedef'd as paramT and attribute of Problem */

class DummyEnergyModel: public EnergyModel {


};

#ifdef USE_VIENNA_ENERGY_MODEL

class ViennaEnergyModel: public EnergyModel {

public:

    /* constructor - initialize the energy parameters here */
    ViennaEnergyModel() : EnergyModel((scaledEnergyParams*) scale_parameters()) {};
    //void update_fold_params() override {
    //    update_fold_params;
    //}
};

#endif

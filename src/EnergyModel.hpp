/*
 * Class for the energy model
 *
 */

extern "C" {
    #include "params.h"
    #include "utils.h"
    #include "energy_const.h"
    #include "fold.h"
}


/*
 * struct for holding RNA scaled energy parameters
 * from ViennaRNA/include/ViennaRNA/params/basic.h vrna_param_s
 */


//// TODO - will need the constants used in these
//typedef struct scaledEnergyParams {
//  int       id;
//  int       stack[NBPAIRS + 1][NBPAIRS + 1];
//  int       hairpin[31];
//  int       bulge[MAXLOOP + 1];
//  int       internal_loop[MAXLOOP + 1];
//  int       mismatchExt[NBPAIRS + 1][5][5];
//  int       mismatchI[NBPAIRS + 1][5][5];
//  int       mismatch1nI[NBPAIRS + 1][5][5];
//  int       mismatch23I[NBPAIRS + 1][5][5];
//  int       mismatchH[NBPAIRS + 1][5][5];
//  int       mismatchM[NBPAIRS + 1][5][5];
//  int       dangle5[NBPAIRS + 1][5];
//  int       dangle3[NBPAIRS + 1][5];
//  int       int11[NBPAIRS + 1][NBPAIRS + 1][5][5];
//  int       int21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
//  int       int22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
//  int       ninio[5];
//  double    lxc;
//  int       MLbase;
//  int       MLintern[NBPAIRS + 1];
//  int       MLclosing;
//  int       TerminalAU;
//  int       DuplexInit;
//  int       Tetraloop_E[200];
//  char      Tetraloops[1401];
//  int       Triloop_E[40];
//  char      Triloops[241];
//  int       Hexaloop_E[40];
//  char      Hexaloops[1801];
//  int       TripleC;
//  int       MultipleCA;
//  int       MultipleCB;
//  int       gquad [VRNA_GQUAD_MAX_STACK_SIZE + 1]
//  [3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
//
//  double    temperature;      /**<  @brief  Temperature used for loop contribution scaling */
//
//  vrna_md_t model_details;    /**<  @brief  Model details to be used in the recursions */
//  char      param_file[256];  /**<  @brief  The filename the parameters were derived from, or empty string if they represent the default */
//
//} scaledEnergyParams;
//
//
///* base class representing an energy model */
//class EnergyModel {
//public:
//    scaledEnergyParams _P;
//
//}
//
//
///* Need something like vrna_param_s type struct for holding temperature scaled energy
// * parameters - typedef'd as paramT and attribute of Problem */
//
//class DummyEnergyModel public EnergyModel {
//
//
//}
//
//
///* Need something like vrna_param_s type struct for holding temperature scaled energy
// * parameters - typedef'd as paramT and attribute of Problem */
//
//class ViennaEnergyModel public EnergyModel {
//
//public:
//
//    /* constructor - initialize the energy parameters here */
//    ViennaEnergyModel() {
//        _P = scale_parameters();
//    }
//
//}

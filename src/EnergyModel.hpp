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
/* vrna_md_t */
struct modelDetails_t {
  double  temperature;                      /**<  @brief  The temperature used to scale the thermodynamic parameters */
  double  betaScale;                        /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
  int     dangles;                          /**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
                                             *
                                             *    If set to 0 no stabilizing energies are assigned to bases adjacent to
                                             *    helices in free ends and multiloops (so called dangling ends). Normally
                                             *    (dangles = 1) dangling end energies are assigned only to unpaired
                                             *    bases and a base cannot participate simultaneously in two dangling ends. In
                                             *    the partition function algorithm vrna_pf() these checks are neglected.
                                             *    To provide comparability between free energy minimization and partition function
                                             *    algorithms, the default setting is 2.
                                             *    This treatment of dangling ends gives more favorable energies to helices
                                             *    directly adjacent to one another, which can be beneficial since such
                                             *    helices often do engage in stabilizing interactions through co-axial
                                             *    stacking.\n
                                             *    If set to 3 co-axial stacking is explicitly included for
                                             *    adjacent helices in multiloops. The option affects only mfe folding
                                             *    and energy evaluation (vrna_mfe() and vrna_eval_structure()), as
                                             *    well as suboptimal folding (vrna_subopt()) via re-evaluation of energies.
                                             *    Co-axial stacking with one intervening mismatch is not considered so far.
                                             *    @note   Some function do not implement all dangle model but only a subset of
                                             *            (0,1,2,3). In particular, partition function algorithms can only handle
                                             *            0 and 2. Read the documentation of the particular recurrences or
                                             *            energy evaluation function for information about the provided dangle
                                             *            model.
                                             */
  int     special_hp;                       /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
  int     noLP;                             /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
  int     noGU;                             /**<  @brief  Do not allow GU pairs */
  int     noGUclosure;                      /**<  @brief  Do not allow loops to be closed by GU pair */
  int     logML;                            /**<  @brief  Use logarithmic scaling for multiloops */
  int     circ;                             /**<  @brief  Assume RNA to be circular instead of linear */
  int     gquad;                            /**<  @brief  Include G-quadruplexes in structure prediction */
  int     uniq_ML;                          /**<  @brief  Flag to ensure unique multi-branch loop decomposition during folding */
  int     energy_set;                       /**<  @brief  Specifies the energy set that defines set of compatible base pairs */
  int     backtrack;                        /**<  @brief  Specifies whether or not secondary structures should be backtraced */
  char    backtrack_type;                   /**<  @brief  Specifies in which matrix to backtrack */
  int     compute_bpp;                      /**<  @brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
  char    nonstandards[64];                 /**<  @brief  contains allowed non standard bases */
  int     max_bp_span;                      /**<  @brief  maximum allowed base pair span */

  int     min_loop_size;                    /**<  @brief  Minimum size of hairpin loops
                                             *    @note The default value for this field is #TURN, however, it may
                                             *    be 0 in cofolding context.
                                             */
  int     window_size;                      /**<  @brief  Size of the sliding window for locally optimal structure prediction */
  int     oldAliEn;                         /**<  @brief  Use old alifold energy model */
  int     ribo;                             /**<  @brief  Use ribosum scoring table in alifold energy model */
  double  cv_fact;                          /**<  @brief  Co-variance scaling factor for consensus structure prediction */
  double  nc_fact;                          /**<  @brief  Scaling factor to weight co-variance contributions of non-canonical pairs */
  double  sfact;                            /**<  @brief  Scaling factor for partition function scaling */
  int     rtype[8];                         /**<  @brief  Reverse base pair type array */
  short   alias[MAXALPHA + 1];              /**<  @brief  alias of an integer nucleotide representation */
  int     pair[MAXALPHA + 1][MAXALPHA + 1]; /**<  @brief  Integer representation of a base pair */
};


/*
 * struct for holding RNA scaled energy parameters
 * from ViennaRNA/include/ViennaRNA/params/basic.h vrna_param_s
 */


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

  modelDetails_t model_details;    /**<  @brief  Model details to be used in the recursions */
  char      param_file[256];  /**<  @brief  The filename the parameters were derived from, or empty string if they represent the default */

} scaledEnergyParams;


/* base class representing an energy model */
class EnergyModel {
protected:
    std::shared_ptr<scaledEnergyParams> energyParams_;  // TODO - make this protected

public:
    /* constructor 
     * inputs: p - pointer to a scaledEnergyParams struct that will be copied
     *             into the energyParams_ attribute */ 
    EnergyModel(scaledEnergyParams* p = nullptr) : energyParams_(p) {};
   
    /* moveEnergyParams - moves the unique pointer to energyParams
     */
    std::shared_ptr<scaledEnergyParams> getEnergyParams () {
        return energyParams_;
    }
    
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
    DummyEnergyModel() : EnergyModel(new scaledEnergyParams) {};

    void repr() override {
        std::cout << "Dummy Energy Model" << std::endl;
    }

    void updateEnergyFoldParams() override {

    }

};

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

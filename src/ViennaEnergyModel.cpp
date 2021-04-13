#include <cstring>
#include <cmath>
#include <memory>

#include "EnergyModel.hpp"
#include "ViennaEnergyModel.hpp"
//extern "C" {
//#include "params.h"
//#include "utils.h"
//}

int ViennaEnergyModel::TermAU(int const &type){
    if (type > 2) {
        return energyParams_->TerminalAU;
    }
    return 0;
}

int ViennaEnergyModel::E_hairpin(int size, int type, int si1, int sj1, const char *string) {
    int energy = (size <= 30) ? energyParams_->hairpin[size] : energyParams_->hairpin[30] + (int)(energyParams_->lxc * log((size) / 30.));
    // fprintf(stderr, "ok\n");
    if (energyParams_->model_details.special_hp) {
        char tl[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, *ts;
        if (size == 4) { /* check for tetraloop bonus */
            strcpy(tl, string);
            if ((ts = strstr(energyParams_->Tetraloops, tl))) {
                return (energyParams_->Tetraloop_E[(ts - energyParams_->Tetraloops) / 7]);
            }
        } else if (size == 6) {
            strcpy(tl, string);
            if ((ts = strstr(energyParams_->Hexaloops, tl))) {
                return (energy = energyParams_->Hexaloop_E[(ts - energyParams_->Hexaloops) / 9]);
            }
        } else if (size == 3) {
            strcpy(tl, string);
            if ((ts = strstr(energyParams_->Triloops, tl))) {
                return (energyParams_->Triloop_E[(ts - energyParams_->Triloops) / 6]);
            }
            return (energy + (type > 2 ? energyParams_->TerminalAU : 0));
        }
    }
    energy += energyParams_->mismatchH[type][si1][sj1];

    return energy;
}

int ViennaEnergyModel::E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1){
    /* compute energy of degree 2 loop (stack bulge or interior) */
    // AMW: trying to reduce the scope of the variable energy surprisingly reduces performance a small amount.
    // Might be an inlining issue.
    int energy = INF;
    int MAX_NINIO = 300;
    int const nl = n1 > n2 ? n1 : n2;
    int const ns = n1 > n2 ? n2 : n1;

    if (nl == 0)
        return energyParams_->stack[type][type_2]; /* stack */

    if (ns == 0) { /* bulge */
        energy = (nl <= MAXLOOP) ? energyParams_->bulge[nl] : (energyParams_->bulge[30] + (int)(energyParams_->lxc * log(nl / 30.)));
        if (nl == 1)
            energy += energyParams_->stack[type][type_2];
        else {
            if (type > 2)
                energy += energyParams_->TerminalAU;
            if (type_2 > 2)
                energy += energyParams_->TerminalAU;
        }
        return energy;
    } else { /* interior loop */
        if (ns == 1) {
            if (nl == 1) /* 1x1 loop */
                return energyParams_->int11[type][type_2][si1][sj1];
            if (nl == 2) { /* 2x1 loop */
                if (n1 == 1)
                    return energyParams_->int21[type][type_2][si1][sq1][sj1];
                else
                    return energyParams_->int21[type_2][type][sq1][si1][sp1];
            } else { /* 1xn loop */
                energy = (nl + 1 <= MAXLOOP) ? (energyParams_->internal_loop[nl + 1])
                                             : (energyParams_->internal_loop[30] + (int)(energyParams_->lxc * log((nl + 1) / 30.)));
                energy += MIN2(MAX_NINIO, (nl - ns) * energyParams_->ninio[2]);
                energy += energyParams_->mismatch1nI[type][si1][sj1] + energyParams_->mismatch1nI[type_2][sq1][sp1];
                return energy;
            }
        } else if (ns == 2) {
            if (nl == 2) { /* 2x2 loop */
                return energyParams_->int22[type][type_2][si1][sp1][sq1][sj1];
            } else if (nl == 3) { /* 2x3 loop */
                energy = energyParams_->internal_loop[5] + energyParams_->ninio[2];
                energy += energyParams_->mismatch23I[type][si1][sj1] + energyParams_->mismatch23I[type_2][sq1][sp1];
                return energy;
            }
        }
        { /* generic interior loop (no else here!)*/
            energy = (n1 + n2 <= MAXLOOP) ? (energyParams_->internal_loop[n1 + n2])
                                          : (energyParams_->internal_loop[30] + (int)(energyParams_->lxc * log((n1 + n2) / 30.)));
            energy += MIN2(MAX_NINIO, (nl - ns) * energyParams_->ninio[2]);

            energy += energyParams_->mismatchI[type][si1][sj1] + energyParams_->mismatchI[type_2][sq1][sp1];
            return energy;
        }
    }
}


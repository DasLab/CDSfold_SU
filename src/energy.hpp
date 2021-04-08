#pragma once

#include <cstring>
#include <cmath>
#include <memory>

#include "EnergyModel.hpp"
//extern "C" {
//#include "params.h"
//#include "utils.h"
//}

inline auto TermAU(int const &type, std::shared_ptr<scaledEnergyParams> const & P) -> int {
    if (type > 2) {
        return P->TerminalAU;
    }
    return 0;
}

inline auto E_hairpin(int size, int type, int si1, int sj1, const char *string, std::shared_ptr<scaledEnergyParams> const & P) -> int {
    int energy = (size <= 30) ? P->hairpin[size] : P->hairpin[30] + (int)(P->lxc * log((size) / 30.));
    // fprintf(stderr, "ok\n");
    if (P->model_details.special_hp) {
        char tl[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, *ts;
        if (size == 4) { /* check for tetraloop bonus */
            strcpy(tl, string);
            if ((ts = strstr(P->Tetraloops, tl))) {
                return (P->Tetraloop_E[(ts - P->Tetraloops) / 7]);
            }
        } else if (size == 6) {
            strcpy(tl, string);
            if ((ts = strstr(P->Hexaloops, tl))) {
                return (energy = P->Hexaloop_E[(ts - P->Hexaloops) / 9]);
            }
        } else if (size == 3) {
            strcpy(tl, string);
            if ((ts = strstr(P->Triloops, tl))) {
                return (P->Triloop_E[(ts - P->Triloops) / 6]);
            }
            return (energy + (type > 2 ? P->TerminalAU : 0));
        }
    }
    energy += P->mismatchH[type][si1][sj1];

    return energy;
}

inline auto E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, std::shared_ptr<scaledEnergyParams> const & P) -> int {
    /* compute energy of degree 2 loop (stack bulge or interior) */
    // AMW: trying to reduce the scope of the variable energy surprisingly reduces performance a small amount.
    // Might be an inlining issue.
    int energy = INF;
    int MAX_NINIO = 300;
    int const nl = n1 > n2 ? n1 : n2;
    int const ns = n1 > n2 ? n2 : n1;

    if (nl == 0)
        return P->stack[type][type_2]; /* stack */

    if (ns == 0) { /* bulge */
        energy = (nl <= MAXLOOP) ? P->bulge[nl] : (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
        if (nl == 1)
            energy += P->stack[type][type_2];
        else {
            if (type > 2)
                energy += P->TerminalAU;
            if (type_2 > 2)
                energy += P->TerminalAU;
        }
        return energy;
    } else { /* interior loop */
        if (ns == 1) {
            if (nl == 1) /* 1x1 loop */
                return P->int11[type][type_2][si1][sj1];
            if (nl == 2) { /* 2x1 loop */
                if (n1 == 1)
                    return P->int21[type][type_2][si1][sq1][sj1];
                else
                    return P->int21[type_2][type][sq1][si1][sp1];
            } else { /* 1xn loop */
                energy = (nl + 1 <= MAXLOOP) ? (P->internal_loop[nl + 1])
                                             : (P->internal_loop[30] + (int)(P->lxc * log((nl + 1) / 30.)));
                energy += MIN2(MAX_NINIO, (nl - ns) * P->ninio[2]);
                energy += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
                return energy;
            }
        } else if (ns == 2) {
            if (nl == 2) { /* 2x2 loop */
                return P->int22[type][type_2][si1][sp1][sq1][sj1];
            } else if (nl == 3) { /* 2x3 loop */
                energy = P->internal_loop[5] + P->ninio[2];
                energy += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
                return energy;
            }
        }
        { /* generic interior loop (no else here!)*/
            energy = (n1 + n2 <= MAXLOOP) ? (P->internal_loop[n1 + n2])
                                          : (P->internal_loop[30] + (int)(P->lxc * log((n1 + n2) / 30.)));
            energy += MIN2(MAX_NINIO, (nl - ns) * P->ninio[2]);

            energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
            return energy;
        }
    }
}


/*
 * JitteredViennaEnergyModel.hpp
 * 
 * Version of ViennaEnergyModel where the energies have some variance
 */

#include "ViennaEnergyModel.hpp"
#pragma once

class JitteredViennaEnergyModel: public ViennaEnergyModel {

private:
    float jitterRange_;
    std::mt19937 randGenerator_;
    std::uniform_real_distribution<> distribution_;

public:
    /* constructor */
    JitteredViennaEnergyModel(float temp, float range, bool fixedSeed): ViennaEnergyModel(temp) {
        jitterRange_ = range;
        std::mt19937 *genp;

        if (fixedSeed) {
            randGenerator_ = std::mt19937(0); 
        }
        else {
            std::random_device r;
            randGenerator_ = std::mt19937(0); 
        }
        distribution_ = std::uniform_real_distribution<>(1 - range, 1 + range);
    }

    /* repr function */

    void repr() override {
        std::cout << "Jittered Vienna Energy Model" << std::endl;
    }

    float generateJitter() {
        return distribution_(randGenerator_); 
    }

    int TermAU(int const &type) {
        return ViennaEnergyModel::TermAU(type);
    }

    int E_hairpin(int size, int type, int si1, int sj1, const char *string) {
        return ViennaEnergyModel::E_hairpin(size, type, si1, sj1, string);
    }

    int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1) {
        return ViennaEnergyModel::E_intloop(n1, n2, type, type_2, si1, sj1, sp1, sq1); 
    }

    /* accessor functions */
    int getTerminalAU() override { return energyParams_->TerminalAU;}
    int getMLbase() override { return energyParams_->MLbase;}
    int getMLintern(unsigned int i) override { return energyParams_->MLintern[i];}
    int getMLclosing() override { return energyParams_->MLclosing;}

};

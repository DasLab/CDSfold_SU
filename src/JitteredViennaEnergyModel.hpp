/*
 * JitteredViennaEnergyModel.hpp
 * 
 * Version of ViennaEnergyModel where the energies have some variance
 */

#pragma once

#include "ViennaEnergyModel.hpp"

class JitteredViennaEnergyModel: public ViennaEnergyModel {

private:
    float jitterRange_;
    std::random_device r_;
    std::mt19937 randGenerator_;
    std::uniform_real_distribution<> distribution_;
    std::map<int, int> jitterCache_;

public:
    /* constructor */
    JitteredViennaEnergyModel(float temp, float range, bool fixedSeed): 
        ViennaEnergyModel(temp),
        jitterRange_(range),
        randGenerator_(fixedSeed ? 0 : r_()),
        distribution_(1 - range, 1 + range) {}

    /* repr function */
    void repr() override {
        std::cout << "Jittered Vienna Energy Model" << std::endl;
    }

    /* jitters an integer value by multiplying it with a random, uniformly 
     * distributed float centered around one */
    int generateJitter(int value) {
        auto cachedVal = jitterCache_.find(value); 
        if (cachedVal != jitterCache_.end()) {
            return cachedVal->second;
        }

        /* value hasn't been jittered before - cache it so that the same 
         * energy parameter always has the same value; backtrace fails if energy 
         * params change */
        float jitter = distribution_(randGenerator_);
        int jitteredVal = jitter * value;
        jitterCache_[value] = jitteredVal;
        
        return jitteredVal;
    }

    int TermAU(int const &type) {
        return generateJitter(ViennaEnergyModel::TermAU(type));
    }

    int E_hairpin(int size, int type, int si1, int sj1, const char *string) {
        return generateJitter(ViennaEnergyModel::E_hairpin(size, type, si1, sj1, string));
    }

    int E_intloop(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1) {
        return generateJitter(ViennaEnergyModel::E_intloop(n1, n2, type, type_2, si1, sj1, sp1, sq1)); 
    }

    /* accessor functions */
    int getTerminalAU() override { return generateJitter(energyParams_->TerminalAU);}
    int getMLbase() override { return generateJitter(energyParams_->MLbase);}
    int getMLintern(unsigned int i) override { return generateJitter(energyParams_->MLintern[i]);}
    int getMLclosing() override { return generateJitter(energyParams_->MLclosing);}

};

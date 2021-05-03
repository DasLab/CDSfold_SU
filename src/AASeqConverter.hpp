/*
 * AASeqConverter.hpp
 *
 *  Created on: 2014/09/11
 *      Author: kamegai
 */

#include "FileManager.hpp"
#include "Util.hpp"
#include "codon.hpp"
#include <algorithm>
#include <ctime>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#pragma once

using namespace std;

enum class Alphabet {
    BASE_ORIGINAL = 0,
    BASE_EXTENDED = 1
};

class AASeqConverter {
  public:
    AASeqConverter() {
        getBaseEnergy();
    }

    ~AASeqConverter() = default;

    auto getBaseEnergy() -> map<string, int> {
        if (baseEnergy.empty()) {
            FileManager fm;
            fm.loadEnergyFile(baseEnergy);
        }
        return baseEnergy;
    }

    auto countNeighborTwoBase(const string &aaseq, const string &exceptedCodons) -> vector<vector<int>>;

    auto countEveryOtherTwoBase(const string &aaseq, const string &exceptedCodons) -> vector<vector<int>>;

    auto getBases(const string &aminoAcid, const string &exceptedCodons, Alphabet const & flag) 
        -> vector<vector<vector<string>>>;

    auto calcQueryBaseEnergy(const string &aminoAcid, string exceptedCodons, Alphabet const & flag)
        -> vector<vector<vector<vector<pair<int, string>>>>>;


  private:
    codon codonTable;
    static const map<string, int> baseNumberMap; 
    static const map<string, int> pairNumberMap;
    map<string, int> baseEnergy;

    auto getBaseNumber(const string &base) -> int {
        auto itr = baseNumberMap.find(base);
        return itr->second;
    }

    auto getPairNumber(const string &twoBases) -> int {
        auto itr = pairNumberMap.find(twoBases);
        return itr->second;
    }

    void initBasePairMap(int size, vector<vector<int>> &map);

    void setBasePairMap(int position, string basePair, vector<vector<int>> &map);

    void initResultBaseVector(int inputBaseLength, int getBaseLength, vector<vector<vector<string>>> &result);

    void initResultExtendedBaseEnergyVector(int inputBaseLength, int getBaseLength,
                                            vector<vector<vector<vector<pair<int, string>>>>> &result) {
        initResultBaseEnergyVectorBase(inputBaseLength, getBaseLength, result, Alphabet::BASE_EXTENDED);
    }
    
    void initResultOriginalBaseEnergyVector(int inputBaseLength, int getBaseLength,
                                            vector<vector<vector<vector<pair<int, string>>>>> &result) {
        initResultBaseEnergyVectorBase(inputBaseLength, getBaseLength, result, Alphabet::BASE_ORIGINAL);
    }

    void initResultBaseEnergyVectorBase(int inputBaseLength, int getBaseLength,
                                        vector<vector<vector<vector<pair<int, string>>>>> &result, 
                                        Alphabet const & flag);


    void getBasesBase(string aminoAcid, unsigned int maxLength, vector<vector<vector<string>>> &result, 
                      Alphabet const & flag, const string &exceptedCodons);


    void addBasesBase(string aminoAcid, int startElement, unsigned int maxLength,
                      vector<vector<vector<string>>> &result, Alphabet const & flag, 
                      const string &exceptedCodons);

    void addExtendedBases(string aminoAcid, int startElement, unsigned int maxLength,
                          vector<vector<vector<string>>> &result, string exceptedCodons) {
        addBasesBase(std::move(aminoAcid), startElement, maxLength, result, 
                     Alphabet::BASE_EXTENDED, std::move(exceptedCodons));
    }

    void addOriginalBases(string aminoAcid, int startElement, unsigned int maxLength,
                          vector<vector<vector<string>>> &result, string exceptedCodons) {
        addBasesBase(std::move(aminoAcid), startElement, maxLength, result, 
                     Alphabet::BASE_ORIGINAL, std::move(exceptedCodons));
    }

    auto calcEachBaseEnergyBase(const string &aminoAcid, vector<vector<vector<string>>> &bases,
                                Alphabet const & flag) -> vector<vector<vector<vector<pair<int, string>>>>>;

    auto calcEachBaseEnergy(string aminoAcid, vector<vector<vector<string>>> &bases, Alphabet const & flag)
        -> vector<vector<vector<vector<pair<int, string>>>>> {
        return calcEachBaseEnergyBase(std::move(aminoAcid), bases, flag);
    }

    auto getSelectedLengthBases(string aminoAcid, unsigned int maxLength, string exceptedCodons, 
                                Alphabet const & flag)
        -> vector<vector<vector<string>>> {
        vector<vector<vector<string>>> result;
        getBasesBase(std::move(aminoAcid), maxLength, result, flag, std::move(exceptedCodons));
        return result;
    }

    /*
     * For debugging
     */
    void printTime(const string &msg, clock_t start, clock_t end) {
        double margin = (double)(end - start) / CLOCKS_PER_SEC;
        cout << msg << "\t" << margin << endl;
    }
};

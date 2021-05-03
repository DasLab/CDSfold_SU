#include "AASeqConverter.hpp"
    
auto AASeqConverter::countNeighborTwoBase(const string &aaseq, const string &exceptedCodons) -> 
    vector<vector<int>> {
    
    vector<vector<int>> result;

    // Initialize result storage map
    int twoPairSize = aaseq.size() * 3 - 1;
    initBasePairMap(twoPairSize, result);

    // Convert a string of amino acid sequences to a sequence
    int position = 0;
    map<string, int> preBaseMap;
    for (char aa : aaseq) {
        // Get a list of amino acid base candidates
        vector<string> baseList = codonTable.getExtendedCodons(aa, exceptedCodons);

        // Create a base sequence candidate sequence
        map<string, int> nextPreBaseMap;
        for (string const & codon : baseList) {
            for (auto const & elem : preBaseMap) {
                string const & preBase = elem.first;
                string basePair = preBase + codon.substr(0, 1);
                setBasePairMap(position, basePair, result);
            }

            setBasePairMap(position + 1, codon.substr(0, 2), result);
            setBasePairMap(position + 2, codon.substr(1, 2), result);
            nextPreBaseMap.insert(make_pair(codon.substr(2, 1), 1));
        }
        position = position + 3;
        preBaseMap = nextPreBaseMap;
    }
    return result;
}

auto AASeqConverter::countEveryOtherTwoBase(const string &aaseq, const string &exceptedCodons) -> 
    vector<vector<int>> {
    
    vector<vector<int>> result;

    // Initialize result storage map
    int twoPairSize = aaseq.size() * 3 - 2;
    initBasePairMap(twoPairSize, result);

    // Convert a string of amino acid sequences to a sequence
    int position = 0;
    map<string, int> preOneBaseMap;
    map<string, int> preTwoBaseMap;
    for (char aa : aaseq) {
        // Get a list of amino acid base candidates
        vector<string> baseList = codonTable.getExtendedCodons(aa, exceptedCodons);

        // Create a base sequence candidate sequence
        map<string, int> nextOneBaseMap;
        map<string, int> nextTwoBaseMap;
        for (string const & codon : baseList) {
            if (!preOneBaseMap.empty()) {
                for (auto const & elem : preTwoBaseMap) {
                    string const & preBase = elem.first;
                    string basePair = preBase + codon.substr(0, 1);
                    setBasePairMap(position - 1, basePair, result);
                }
                for (auto const & elem : preOneBaseMap) {
                    string const &  preBase = elem.first;
                    string basePair = preBase + codon.substr(1, 1);
                    setBasePairMap(position, basePair, result);
                }
            }
            string twoBase = codon.substr(0, 1) + codon.substr(2, 1);
            setBasePairMap(position + 1, twoBase, result);
            nextTwoBaseMap.insert(make_pair(codon.substr(1, 1), 1));
            nextOneBaseMap.insert(make_pair(codon.substr(2, 1), 1));
        }
        position = position + 3;
        preTwoBaseMap = nextTwoBaseMap;
        preOneBaseMap = nextOneBaseMap;
    }

    return result;
}


auto AASeqConverter::getBases(const string &aminoAcid, const string &exceptedCodons, Alphabet const & flag) -> vector<vector<vector<string>>> {
    // Obtain the partial base sequence of 0 to (n-8) elements
    unsigned int maxLength = 8;
    vector<vector<vector<string>>> bases = getSelectedLengthBases(aminoAcid, maxLength, exceptedCodons, flag);

    // Obtain the partial base sequence after (n-8 + 1)
    int preMaxLength = maxLength;
    for (int i = maxLength - 1; i > 0; --i) {
        maxLength = i;
        int startElement = aminoAcid.size() * 3 - preMaxLength + 1;
        // Actually does the same thing whether Extended or Original!!
        addExtendedBases(aminoAcid, startElement, maxLength, bases, exceptedCodons);
        preMaxLength = maxLength;
    }

    return bases;
}

/*
 * Obtain the free energy of the partial base sequence site of the amino acid sequence
 *
 * @param　aminoAcid　Amino acid sequence
 * @return　Free energy vector in each base sequence
 * 　Element 1D: Start position (number of bases from 0 to amino acid sequence, actually use element number 1 or later)
 * 　Element 2D: Partial array length (0-8, actually element numbers are only 5,6,8 used)
 * 　Element 3D: Start position base number (0-4, actually use element number 1 or later)
 * 　Element 4D: End position base number (0-4, actually element number 1 or later is used)
 * 　　(Base number: A = 1, C = 2, G = 3, U = 4)
 *
 * (Example of use)
 *   string seq = "MLYF";
 * 	 AASeqConverter conv;
 *	 vector<vector<vector<vector<pair<int, string> > > > > result
 *	   = conv.calcQueryBaseEnergy(seq);
 *	 vector<int> baseLengths;
 *	 baseLengths.push_back(5);
 *	 baseLengths.push_back(6);
 *	 baseLengths.push_back(8);
 *	 for (unsigned int start = 1; start <= result.size(); start++) {
 *	 	for (unsigned int size = 0; size < baseLengths.size(); size++) {
 *	 		int baseSize = baseLengths[size];
 *	 		for (int sBase = 1; sBase <= 4; sBase++) {
 *	 			for (int eBase = 1; eBase <= 4; eBase++) {
 *	 				int energy = result[start][baseSize][sBase][eBase].first;
 *	 				string seq = result[start][baseSize][sBase][eBase].second;
 *	 				cout << start << "\t" << baseSize << "\t" << sBase << "\t" << eBase
 *	 						<< "\t" << energy << "\t" << seq << endl;
 *	 			}
 *	 		}
 *	 	}
 *	 }
 */

auto AASeqConverter::calcQueryBaseEnergy(const string &aminoAcid, string exceptedCodons, Alphabet const & flag)
    -> vector<vector<vector<vector<pair<int, string>>>>> {
    // Get partial base sequence
    vector<vector<vector<string>>> bases = getBases(aminoAcid, std::move(exceptedCodons), flag);

    // Get the minimum energy of each part
    vector<vector<vector<vector<pair<int, string>>>>> result = calcEachBaseEnergy(aminoAcid, bases, flag);

    return result;
}


void AASeqConverter::initBasePairMap(int size, vector<vector<int>> &map) {
    // Create an array of columns to be assigned to each row (actually used after element 1) 
    vector<int> rows(size+1, 0);

    // Assign a column to each row (actually used after element 1)
    for (unsigned int i = 0; i <= pairNumberMap.size(); i++) {
        map.push_back(rows);
    }
}

void AASeqConverter::setBasePairMap(int position, string basePair, vector<vector<int>> &map) {
    int pair = getPairNumber(std::move(basePair));
    map[pair][position] = 1;
}

void AASeqConverter::initResultBaseVector(int inputBaseLength, int getBaseLength, vector<vector<vector<string>>> &result) {
    vector<vector<string>> result_length;
    vector<string> result_base;

    for (int i = 0; i <= getBaseLength; i++) {
        result_length.push_back(result_base);
    }

    for (int i = 0; i <= inputBaseLength; i++) {
        result.push_back(result_length);
    }
}


void AASeqConverter::initResultBaseEnergyVectorBase(int inputBaseLength, int getBaseLength,
                                    vector<vector<vector<vector<pair<int, string>>>>> &result, Alphabet const & flag) {
    pair<int, string> baseEnergy;
    baseEnergy.first = numeric_limits<int>::max();
    baseEnergy.second = "";
    // Number of bases
    int baseNumber = 4;
    if (flag == Alphabet::BASE_ORIGINAL) {
        baseNumber = 4;
    } else if (flag == Alphabet::BASE_EXTENDED) {
        baseNumber = 8;
    }
    vector<pair<int, string>> resEndBase;
    vector<vector<pair<int, string>>> resStartBase;
    vector<vector<vector<pair<int, string>>>> resLength;

    // Add an appropriate number of baseEnergy to resEndBase (actually use element 1 or later)
    resEndBase.assign(baseNumber + 1, baseEnergy);

    // Add an appropriate number of resEndBase to resStartBase (actually use element 1 and later)
    resStartBase.assign(baseNumber + 1, resEndBase);

    // Add an appropriate number of resStartBase to resLength (actually use element 1 and above
    resLength.assign(getBaseLength + 1, resStartBase);

    // Add an appropriate number of resLength to result (actually use element 1 and later)
    result.assign(inputBaseLength + 1, resLength);
}

void AASeqConverter::getBasesBase(string aminoAcid, unsigned int maxLength, vector<vector<vector<string>>> &result, Alphabet const & flag,
                  const string &exceptedCodons) {
    int baseLength = aminoAcid.size() * 3;

    // Initialize the result storage vector
    initResultBaseVector(baseLength, maxLength, result);

    // Create a base sequence with a length of maxLength
    for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
        vector<string> bases; // Stores the created base sequence
        unsigned int startPosition = i + 1;
        for (unsigned int j = 0; j < aminoAcid.size(); j++) {
            if (startPosition > (j + 1) * 3) {
                // If it is before the start position, skip the process
                continue;
            } else {
                // Obtain the base sequence linked to the base sequence being created from the codon
                map<string, int> baseMap;
                int nowPosition = (j + 1) * 3 - 2;
                int addLength;
                if (bases.empty()) {
                    addLength = maxLength;
                } else {
                    addLength = maxLength - bases[0].length();
                }

                vector<string> codons;
                if (flag == Alphabet::BASE_ORIGINAL) {
                    codons = codonTable.getCodons(aminoAcid.at(j), exceptedCodons);
                } else if (flag == Alphabet::BASE_EXTENDED) {
                    codons = codonTable.getExtendedCodons(aminoAcid.at(j), exceptedCodons);
                } else {
                    cout << "Error: getBasesBase () flag is an invalid value:" << static_cast<int>(flag) << endl;
                    exit(1);
                }

                for (auto codon : codons) {
                    int addNumber = 0;
                    string addingBase;
                    for (unsigned int l = 0; l < codon.size(); l++) {
                        if (nowPosition + l >= startPosition) {
                            // Skip processing if maxLength is exceeded 
                            if (addLength - addNumber > 0) {
                                addingBase += codon.at(l);
                                addNumber++;
                            }
                        }
                    }

                    baseMap.insert(make_pair(addingBase, 1));
                }

                // Concatenate the base sequence obtained from the codon to the base sequence being created
                vector<string> newBases;
                map<string, int>::iterator itr;
                for (itr = baseMap.begin(); itr != baseMap.end(); itr++) {
                    string base = itr->first;

                    if (bases.empty()) {
                        newBases.push_back(base);
                    } else {
                        for (auto newBase : bases) {
                            newBase += base;
                            newBases.push_back(newBase);
                        }
                    }
                }

                // Replace the base sequence with the new extended result
                if (!newBases.empty()) {
                    bases = newBases;
                }

                // Exit the loop when the partial array length reaches maxLength
                if (bases[0].length() == maxLength) {
                    break;
                }
            }
        }

        // Set the created base sequence in the result storage vector
        for (auto &base : bases) {
            result[startPosition][maxLength].push_back(base);
        }
    }

    // Create a subarray below maxLength
    for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
        int startPosition = i + 1;
        for (int j = maxLength; j > 1; j--) {
            map<string, int> baseMap;
            for (string const & base : result[startPosition][j]) {
                baseMap.insert(make_pair(base.substr(0, j - 1), 1));
            }

            for (auto const & elem : baseMap) {
                string const & base = elem.first;
                result[startPosition][j - 1].push_back(base);
            }
        }
    }
}

void AASeqConverter::addBasesBase(string aminoAcid, int startElement, unsigned int maxLength,
                  vector<vector<vector<string>>> &result, Alphabet const & flag, const string &exceptedCodons) {
    int baseLength = aminoAcid.size() * 3;

    // Create a base sequence with a length of maxLength
    vector<vector<vector<string>>> nowResult;
    initResultBaseVector(baseLength, maxLength, nowResult);
    for (unsigned int i = startElement; i <= baseLength - maxLength; i++) {
        vector<string> bases; // Stores the created base sequence
        unsigned int startPosition = i + 1;
        for (unsigned int j = 0; j < aminoAcid.size(); j++) {
            if (startPosition > (j + 1) * 3) {
                // If it is before the start position, skip the process
                continue;
            } else {
                // Obtain the base sequence linked to the base sequence being created from the codon
                map<string, int> baseMap;
                int nowPosition = (j + 1) * 3 - 2;
                int addLength;
                if (bases.empty()) {
                    addLength = maxLength;
                } else {
                    addLength = maxLength - bases[0].length();
                }

                vector<string> codons;
                if (flag == Alphabet::BASE_ORIGINAL) {
                    codons = codonTable.getCodons(aminoAcid.at(j), exceptedCodons);
                } else if (flag == Alphabet::BASE_EXTENDED) {
                    codons = codonTable.getExtendedCodons(aminoAcid.at(j), exceptedCodons);
                } else {
                    cout << "Error: The flag in addBasesBase () is an invalid value:" << static_cast<int>(flag) << endl;
                    exit(1);
                }

                for (auto codon : codons) {
                    int addNumber = 0;
                    string addingBase;
                    for (unsigned int l = 0; l < codon.size(); l++) {
                        if (nowPosition + l >= startPosition) {
                            // Skip processing if maxLength is exceeded
                            if (addLength - addNumber > 0) {
                                addingBase += codon.at(l);
                                addNumber++;
                            }
                        }
                    }

                    baseMap.insert(make_pair(addingBase, 1));
                }

                // Concatenate the base sequence obtained from the codon to the base sequence being created
                vector<string> newBases;
                for (auto const & elem : baseMap) {
                    string const & base = elem.first;

                    if (bases.empty()) {
                        newBases.push_back(base);
                    } else {
                        for (auto newBase : bases) {
                            newBase += base;
                            newBases.push_back(newBase);
                        }
                    }
                }

                // Replace the base sequence with the new extended result
                if (!newBases.empty()) {
                    bases = newBases;
                }

                // Exit the loop when the partial array length reaches maxLength
                if (bases[0].length() == maxLength) {
                    break;
                }
            }
        }

        // Set the created base sequence in the result storage vector
        for (auto const & base : bases) {
            nowResult[startPosition][maxLength].push_back(base);
        }
    }

    // Create a subarray below maxLength
    for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
        int startPosition = i + 1;
        for (int j = maxLength; j > 1; j--) {
            map<string, int> baseMap;
            for (string const & base : nowResult[startPosition][j]) {
                baseMap.insert(make_pair(base.substr(0, j - 1), 1));
            }

            for (auto const & elem : baseMap) {
                string const & base = elem.first;
                nowResult[startPosition][j - 1].push_back(base);
            }
        }
    }

    // Set the acquired contents in the result vector
    for (unsigned int i = 1; i < nowResult.size(); i++) {
        for (unsigned int j = 1; j < nowResult[i].size(); j++) {
            for (unsigned int k = 0; k < nowResult[i][j].size(); k++) {
                string base = nowResult[i][j][k];
                result[i][j].push_back(base);
            }
        }
    }
}


auto AASeqConverter::calcEachBaseEnergyBase(const string &aminoAcid, vector<vector<vector<string>>> &bases, Alphabet const & flag)
    -> vector<vector<vector<vector<pair<int, string>>>>> {
    // Initialize the result storage vector
    int baseLength = aminoAcid.size() * 3;
    int maxLength = 8;
    vector<vector<vector<vector<pair<int, string>>>>> result;
    initResultExtendedBaseEnergyVector(baseLength, maxLength, result);

    // Obtain the free energy of the partial array at each position (array length = 5, 6, 8)
    vector<int> baseLengths;
    baseLengths.push_back(5);
    baseLengths.push_back(6);
    baseLengths.push_back(8);
    for (int baseSize : baseLengths) {
        for (int i = 0; i < baseLength; i++) {
            int position = i + 1;
            for (auto const & base : bases[position][baseSize]) {

                // Search for free energy from arrays and energy maps
                string originalBase = base;
                if (flag == Alphabet::BASE_EXTENDED) {
                    // When using extended bases, V and W are converted to U, and X and Y are converted to G.
                    Util::baseReplace(originalBase, "V", "U");
                    Util::baseReplace(originalBase, "W", "U");
                    Util::baseReplace(originalBase, "X", "G");
                    Util::baseReplace(originalBase, "Y", "G");
                }

                // Energy acquisition (skip subsequent processing if not set)
                auto itr = baseEnergy.find(originalBase);
                if (itr == baseEnergy.end()) {
                    continue;
                }
                int energy = itr->second;

                // Store in result vector
                int startBase = getBaseNumber(base.substr(0, 1));
                int endBase = getBaseNumber(base.substr(base.length() - 1, 1));
                int resultEnergy = result[position][baseSize][startBase][endBase].first;
                string resultBase = result[position][baseSize][startBase][endBase].second;

                if (resultBase == "") {
                    result[position][baseSize][startBase][endBase].first = energy;
                    result[position][baseSize][startBase][endBase].second = base;
                } else {
                    if (energy < resultEnergy) {
                        result[position][baseSize][startBase][endBase].first = energy;
                        result[position][baseSize][startBase][endBase].second = base;
                    }
                }
            }
        }
    }
    return result;
}

const map<string, int> AASeqConverter::baseNumberMap {
    {"A", 1},
    {"C", 2},
    {"G", 3},
    {"U", 4},
    {"V", 5},
    {"W", 6},
    {"X", 7},
    {"Y", 8}
};

const map<string, int> AASeqConverter::pairNumberMap{
    {"AA", 1},
    {"AC", 2},
    {"AG", 3},
    {"AU", 4},
    {"AV", 5},
    {"AW", 6},
    {"AX", 7},
    {"AY", 8},
    {"CA", 9},
    {"CC", 10},
    {"CG", 11},
    {"CU", 12},
    {"CV", 13},
    {"CW", 14},
    {"CX", 15},
    {"CY", 16},
    {"GA", 17},
    {"GC", 18},
    {"GG", 19},
    {"GU", 20},
    {"GV", 21},
    {"GW", 22},
    {"GX", 23},
    {"GY", 24},
    {"UA", 25},
    {"UC", 26},
    {"UG", 27},
    {"UU", 28},
    {"UV", 29},
    {"UW", 30},
    {"UX", 31},
    {"UY", 32},
    {"VA", 33},
    {"VC", 34},
    {"VG", 35},
    {"VU", 36},
    {"VV", 37},
    {"VW", 38},
    {"VX", 39},
    {"VY", 40},
    {"WA", 41},
    {"WC", 42},
    {"WG", 43},
    {"WU", 44},
    {"WV", 45},
    {"WW", 46},
    {"WX", 47},
    {"WY", 48},
    {"XA", 49},
    {"XC", 50},
    {"XG", 51},
    {"XU", 52},
    {"XV", 53},
    {"XW", 54},
    {"XX", 55},
    {"XY", 56},
    {"YA", 57},
    {"YC", 58},
    {"YG", 59},
    {"YU", 60},
    {"YV", 61},
    {"YW", 62},
    {"YX", 63},
    {"YY", 64}
};

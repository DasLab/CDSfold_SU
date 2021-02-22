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

    auto countNeighborTwoBase(const string &aaseq, const string &exceptedCodons) -> vector<vector<int>> {
        vector<vector<int>> result;

        // 結果格納マップを初期化
        int twoPairSize = aaseq.size() * 3 - 1;
        initBasePairMap(twoPairSize, result);

        // アミノ酸配列の文字列を配列に変換
        int position = 0;
        map<string, int> preBaseMap;
        for (char aa : aaseq) {
            // アミノ酸の塩基候補リストを取得
            vector<string> baseList = codonTable.getExtendedCodons(aa, exceptedCodons);

            // 塩基配列の候補配列を作成
            map<string, int> nextPreBaseMap;
            vector<string>::iterator itr;
            itr = baseList.begin();
            while (itr != baseList.end()) {
                string codon = *itr;

                if (preBaseMap.empty()) {

                } else {
                    map<string, int>::iterator mapItr;
                    for (mapItr = preBaseMap.begin(); mapItr != preBaseMap.end(); mapItr++) {
                        string preBase = mapItr->first;
                        string basePair = preBase + codon.substr(0, 1);
                        setBasePairMap(position, basePair, result);
                    }
                }

                setBasePairMap(position + 1, codon.substr(0, 2), result);
                setBasePairMap(position + 2, codon.substr(1, 2), result);
                nextPreBaseMap.insert(make_pair(codon.substr(2, 1), 1));
                itr++;
            }
            position = position + 3;
            preBaseMap = nextPreBaseMap;
        }

        return result;
    }

    auto countEveryOtherTwoBase(const string &aaseq, const string &exceptedCodons) -> vector<vector<int>> {
        vector<vector<int>> result;

        // 結果格納マップを初期化
        int twoPairSize = aaseq.size() * 3 - 2;
        initBasePairMap(twoPairSize, result);

        // アミノ酸配列の文字列を配列に変換
        int position = 0;
        map<string, int> preOneBaseMap;
        map<string, int> preTwoBaseMap;
        for (char aa : aaseq) {
            // アミノ酸の塩基候補リストを取得
            vector<string> baseList = codonTable.getExtendedCodons(aa, exceptedCodons);

            // 塩基配列の候補配列を作成
            map<string, int> nextOneBaseMap;
            map<string, int> nextTwoBaseMap;
            vector<string>::iterator itr;
            itr = baseList.begin();
            while (itr != baseList.end()) {
                string codon = *itr;

                if (preOneBaseMap.empty()) {

                } else {
                    map<string, int>::iterator mapItr;
                    for (mapItr = preTwoBaseMap.begin(); mapItr != preTwoBaseMap.end(); mapItr++) {
                        string preBase = mapItr->first;
                        string basePair = preBase + codon.substr(0, 1);
                        setBasePairMap(position - 1, basePair, result);
                    }
                    for (mapItr = preOneBaseMap.begin(); mapItr != preOneBaseMap.end(); mapItr++) {
                        string preBase = mapItr->first;
                        string basePair = preBase + codon.substr(1, 1);
                        setBasePairMap(position, basePair, result);
                    }
                }
                string twoBase = codon.substr(0, 1) + codon.substr(2, 1);
                setBasePairMap(position + 1, twoBase, result);
                nextTwoBaseMap.insert(make_pair(codon.substr(1, 1), 1));
                nextOneBaseMap.insert(make_pair(codon.substr(2, 1), 1));
                itr++;
            }
            position = position + 3;
            preTwoBaseMap = nextTwoBaseMap;
            preOneBaseMap = nextOneBaseMap;
        }

        return result;
    }

    auto getExtendedBases(const string &aminoAcid, const string &exceptedCodons) -> vector<vector<vector<string>>> {
        // 0～(n-8)要素の部分塩基配列を取得
        unsigned int maxLength = 8;
        vector<vector<vector<string>>> bases = getSelectedLengthExtendedBases(aminoAcid, maxLength, exceptedCodons);

        // (n-8+1)以降の部分塩基配列を取得
        int preMaxLength = maxLength;
        for (int i = maxLength - 1; i > 0; --i) {
            maxLength = i;
            int startElement = aminoAcid.size() * 3 - preMaxLength + 1;
            addExtendedBases(aminoAcid, startElement, maxLength, bases, exceptedCodons);
            preMaxLength = maxLength;
        }

        return bases;
    }

    /*
     * アミノ酸配列の部分塩基配列部位の自由エネルギーを取得する
     *
     * @param　aminoAcid　アミノ酸配列
     * @return　各塩基配列における自由エネルギーのベクター
     * 　要素1次元：開始位置（0～アミノ酸配列の塩基数、実際には要素番号1以降を使用）
     * 　要素2次元：部分配列長（0～8、実際には要素番号は、5,6,8のみ使用）
     * 　要素3次元：開始位置塩基番号（0～8、実際には要素番号1以降を使用）
     * 　要素4次元：終了位置塩基番号（0～8、実際には要素番号1以降を使用）
     * 　　（塩基番号：A=1, C=2, G=3, U=4, V=5, W=6, X=7, Y=8)
     *
     * （使用例）
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
     *	 		for (int sBase = 1; sBase <= 8; sBase++) {
     *	 			for (int eBase = 1; eBase <= 8; eBase++) {
     *	 				int energy = result[start][baseSize][sBase][eBase].first;
     *	 				string seq = result[start][baseSize][sBase][eBase].second;
     *	 				cout << start << "\t" << baseSize << "\t" << sBase << "\t" << eBase
     *	 						<< "\t" << energy << "\t" << seq << endl;
     *	 			}
     *	 		}
     *	 	}
     *	 }
     */
    auto calcQueryExtendedBaseEnergy(const string &aminoAcid, string exceptedCodons)
        -> vector<vector<vector<vector<pair<int, string>>>>> {
        // 部分配列を取得
        vector<vector<vector<string>>> bases = getExtendedBases(aminoAcid, std::move(exceptedCodons));

        // 各部位の最小エネルギーを取得
        vector<vector<vector<vector<pair<int, string>>>>> result = calcEachExtendedBaseEnergy(aminoAcid, bases);

        return result;
    }

    auto getOriginalBases(const string &aminoAcid, const string &exceptedCodons) -> vector<vector<vector<string>>> {
        // 0～(n-8)要素の部分塩基配列を取得
        unsigned int maxLength = 8;
        vector<vector<vector<string>>> bases = getSelectedLengthOriginalBases(aminoAcid, maxLength, exceptedCodons);

        // (n-8+1)以降の部分塩基配列を取得
        int preMaxLength = maxLength;
        for (int i = maxLength - 1; i > 0; --i) {
            maxLength = i;
            int startElement = aminoAcid.size() * 3 - preMaxLength + 1;
            addExtendedBases(aminoAcid, startElement, maxLength, bases, exceptedCodons);
            preMaxLength = maxLength;
        }

        //		// (n-8+1)～(n-6)要素の部分塩基配列を取得
        //		int preMaxLength = 8;
        //		int maxLength = 6;
        //		int startElement = aminoAcid.size() * 3 - preMaxLength + 1;
        //		addOriginalBases(aminoAcid, startElement, maxLength, bases,
        //				exceptedCodons);
        //
        //		// (n-6+1)～(n-5)要素の部分塩基配列を取得
        //		preMaxLength = maxLength;
        //		maxLength = 5;
        //		startElement = aminoAcid.size() * 3 - preMaxLength + 1;
        //		addOriginalBases(aminoAcid, startElement, maxLength, bases,
        //				exceptedCodons);

        return bases;
    }

    /*
     * アミノ酸配列の部分塩基配列部位の自由エネルギーを取得する
     *
     * @param　aminoAcid　アミノ酸配列
     * @return　各塩基配列における自由エネルギーのベクター
     * 　要素1次元：開始位置（0～アミノ酸配列の塩基数、実際には要素番号1以降を使用）
     * 　要素2次元：部分配列長（0～8、実際には要素番号は、5,6,8のみ使用）
     * 　要素3次元：開始位置塩基番号（0～4、実際には要素番号1以降を使用）
     * 　要素4次元：終了位置塩基番号（0～4、実際には要素番号1以降を使用）
     * 　　（塩基番号：A=1, C=2, G=3, U=4)
     *
     * （使用例）
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
    auto calcQueryOriginalBaseEnergy(const string &aminoAcid, string exceptedCodons)
        -> vector<vector<vector<vector<pair<int, string>>>>> {
        // 部分塩基配列を取得
        vector<vector<vector<string>>> bases = getOriginalBases(aminoAcid, std::move(exceptedCodons));

        // 各部位の最小エネルギーを取得
        vector<vector<vector<vector<pair<int, string>>>>> result = calcEachOriginalBaseEnergy(aminoAcid, bases);

        return result;
    }

  private:
    codon codonTable;
    static const map<string, int> baseNumberMap; 
    static const map<string, int> pairNumberMap;
    map<string, int> baseEnergy;

    const static int BASE_ORIGINAL = 0;
    const static int BASE_EXTENDED = 1;

    auto getBaseNumber(const string &base) -> int {
        auto itr = baseNumberMap.find(base);
        return itr->second;
    }

    auto getPairNumber(const string &twoBases) -> int {
        auto itr = pairNumberMap.find(twoBases);
        return itr->second;
    }

    void initBasePairMap(int size, vector<vector<int>> &map) {
        // 各行に代入する列の配列を作成（実際に使用するのは、要素1以降）
        vector<int> rows;
        for (int i = 0; i <= size; i++) {
            rows.push_back(0);
        }

        // 各行に列を代入（実際に使用するのは、要素1以降）
        for (unsigned int i = 0; i <= pairNumberMap.size(); i++) {
            map.push_back(rows);
        }
    }

    void setBasePairMap(int position, string basePair, vector<vector<int>> &map) {
        int pair = getPairNumber(std::move(basePair));
        map[pair][position] = 1;
    }

    void initResultBaseVector(int inputBaseLength, int getBaseLength, vector<vector<vector<string>>> &result) {
        vector<vector<string>> result_length;
        vector<string> result_base;

        for (int i = 0; i <= getBaseLength; i++) {
            result_length.push_back(result_base);
        }

        for (int i = 0; i <= inputBaseLength; i++) {
            result.push_back(result_length);
        }
    }

    void initResultExtendedBaseEnergyVector(int inputBaseLength, int getBaseLength,
                                            vector<vector<vector<vector<pair<int, string>>>>> &result) {
        initResultBaseEnergyVectorBase(inputBaseLength, getBaseLength, result, BASE_EXTENDED);
    }

    void initResultOriginalBaseEnergyVector(int inputBaseLength, int getBaseLength,
                                            vector<vector<vector<vector<pair<int, string>>>>> &result) {
        initResultBaseEnergyVectorBase(inputBaseLength, getBaseLength, result, BASE_ORIGINAL);
    }

    void initResultBaseEnergyVectorBase(int inputBaseLength, int getBaseLength,
                                        vector<vector<vector<vector<pair<int, string>>>>> &result, int flag) {
        pair<int, string> baseEnergy;
        baseEnergy.first = numeric_limits<int>::max();
        baseEnergy.second = "";
        // 塩基数
        int baseNumber;
        if (flag == BASE_ORIGINAL) {
            baseNumber = 4;
        } else if (flag == BASE_EXTENDED) {
            baseNumber = 8;
        }
        vector<pair<int, string>> resEndBase;
        vector<vector<pair<int, string>>> resStartBase;
        vector<vector<vector<pair<int, string>>>> resLength;

        // resEndBaseに適切な数のbaseEnergyを追加（実際に使用するのは、要素1以降）
        resEndBase.assign(baseNumber + 1, baseEnergy);

        // resStartBaseに適切な数のresEndBaseを追加（実際に使用するのは、要素1以降）
        resStartBase.assign(baseNumber + 1, resEndBase);

        // resLengthに適切な数のresStartBaseを追加（実際に使用するのは、要素1以降）
        resLength.assign(getBaseLength + 1, resStartBase);

        // resultに適切な数のresLengthを追加（実際に使用するのは、要素1以降）
        result.assign(inputBaseLength + 1, resLength);
    }

    void getBasesBase(string aminoAcid, unsigned int maxLength, vector<vector<vector<string>>> &result, int flag,
                      const string &exceptedCodons) {
        int baseLength = aminoAcid.size() * 3;

        // 結果格納ベクターを初期化
        initResultBaseVector(baseLength, maxLength, result);

        // maxLengthの長さを持つ塩基配列を作成
        for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
            vector<string> bases; // 作成した塩基配列を格納
            unsigned int startPosition = i + 1;
            for (unsigned int j = 0; j < aminoAcid.size(); j++) {
                if (startPosition > (j + 1) * 3) {
                    // 開始位置より手前なら処理をスキップ
                    continue;
                } else {
                    // 作成中の塩基配列に連結する塩基配列をコドンから取得
                    map<string, int> baseMap;
                    int nowPosition = (j + 1) * 3 - 2;
                    int addLength;
                    if (bases.empty()) {
                        addLength = maxLength;
                    } else {
                        addLength = maxLength - bases[0].length();
                    }

                    vector<string> codons;
                    if (flag == BASE_ORIGINAL) {
                        codons = codonTable.getCodons(aminoAcid.at(j), exceptedCodons);
                    } else if (flag == BASE_EXTENDED) {
                        codons = codonTable.getExtendedCodons(aminoAcid.at(j), exceptedCodons);
                    } else {
                        cout << "エラー：getBasesBase()のflagが無効な値です:" << flag << endl;
                        exit(1);
                    }

                    for (auto codon : codons) {
                        int addNumber = 0;
                        string addingBase;
                        for (unsigned int l = 0; l < codon.size(); l++) {
                            if (nowPosition + l >= startPosition) {
                                // maxLengthを超える場合は処理をスキップ
                                if (addLength - addNumber > 0) {
                                    addingBase += codon.at(l);
                                    addNumber++;
                                }
                            }
                        }

                        baseMap.insert(make_pair(addingBase, 1));
                    }

                    // コドンから取得した塩基配列を作成中の塩基配列に連結
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

                    // 塩基配列を伸長した新しい結果に置き換える
                    if (!newBases.empty()) {
                        bases = newBases;
                    }

                    // 部分配列の長さmaxLengthに達したらループを抜ける
                    if (bases[0].length() == maxLength) {
                        break;
                    }
                }
            }

            // 作成した塩基配列を結果格納ベクターにセット
            for (auto &base : bases) {
                result[startPosition][maxLength].push_back(base);
            }
        }

        // maxLength以下の部分配列を作成
        for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
            int startPosition = i + 1;
            for (int j = maxLength; j > 1; j--) {
                map<string, int> baseMap;
                for (unsigned int k = 0; k < result[startPosition][j].size(); k++) {
                    string base = result[startPosition][j][k];
                    baseMap.insert(make_pair(base.substr(0, j - 1), 1));
                }

                map<string, int>::iterator itr;
                for (itr = baseMap.begin(); itr != baseMap.end(); itr++) {
                    string base = itr->first;
                    result[startPosition][j - 1].push_back(base);
                }
            }
        }
    }

    void addBasesBase(string aminoAcid, int startElement, unsigned int maxLength,
                      vector<vector<vector<string>>> &result, int flag, const string &exceptedCodons) {
        int baseLength = aminoAcid.size() * 3;

        // maxLengthの長さを持つ塩基配列を作成
        vector<vector<vector<string>>> nowResult;
        initResultBaseVector(baseLength, maxLength, nowResult);
        for (unsigned int i = startElement; i <= baseLength - maxLength; i++) {
            vector<string> bases; // 作成した塩基配列を格納
            unsigned int startPosition = i + 1;
            for (unsigned int j = 0; j < aminoAcid.size(); j++) {
                if (startPosition > (j + 1) * 3) {
                    // 開始位置より手前なら処理をスキップ
                    continue;
                } else {
                    // 作成中の塩基配列に連結する塩基配列をコドンから取得
                    map<string, int> baseMap;
                    int nowPosition = (j + 1) * 3 - 2;
                    int addLength;
                    if (bases.empty()) {
                        addLength = maxLength;
                    } else {
                        addLength = maxLength - bases[0].length();
                    }

                    vector<string> codons;
                    if (flag == BASE_ORIGINAL) {
                        codons = codonTable.getCodons(aminoAcid.at(j), exceptedCodons);
                    } else if (flag == BASE_EXTENDED) {
                        codons = codonTable.getExtendedCodons(aminoAcid.at(j), exceptedCodons);
                    } else {
                        cout << "エラー：addBasesBase()のflagが無効な値です:" << flag << endl;
                        exit(1);
                    }

                    for (auto codon : codons) {
                        int addNumber = 0;
                        string addingBase;
                        for (unsigned int l = 0; l < codon.size(); l++) {
                            if (nowPosition + l >= startPosition) {
                                // maxLengthを超える場合は処理をスキップ
                                if (addLength - addNumber > 0) {
                                    addingBase += codon.at(l);
                                    addNumber++;
                                }
                            }
                        }

                        baseMap.insert(make_pair(addingBase, 1));
                    }

                    // コドンから取得した塩基配列を作成中の塩基配列に連結
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

                    // 塩基配列を伸長した新しい結果に置き換える
                    if (!newBases.empty()) {
                        bases = newBases;
                    }

                    // 部分配列の長さmaxLengthに達したらループを抜ける
                    if (bases[0].length() == maxLength) {
                        break;
                    }
                }
            }

            // 作成した塩基配列を結果格納ベクターにセット
            for (auto &base : bases) {
                nowResult[startPosition][maxLength].push_back(base);
            }
        }

        // maxLength以下の部分配列を作成
        for (unsigned int i = 0; i <= baseLength - maxLength; i++) {
            int startPosition = i + 1;
            for (int j = maxLength; j > 1; j--) {
                map<string, int> baseMap;
                for (unsigned int k = 0; k < nowResult[startPosition][j].size(); k++) {
                    string base = nowResult[startPosition][j][k];
                    baseMap.insert(make_pair(base.substr(0, j - 1), 1));
                }

                map<string, int>::iterator itr;
                for (itr = baseMap.begin(); itr != baseMap.end(); itr++) {
                    string base = itr->first;
                    nowResult[startPosition][j - 1].push_back(base);
                }
            }
        }

        // 取得した内容を結果ベクターに設定
        for (unsigned int i = 1; i < nowResult.size(); i++) {
            for (unsigned int j = 1; j < nowResult[i].size(); j++) {
                for (unsigned int k = 0; k < nowResult[i][j].size(); k++) {
                    string base = nowResult[i][j][k];
                    result[i][j].push_back(base);
                }
            }
        }
    }

    void addExtendedBases(string aminoAcid, int startElement, unsigned int maxLength,
                          vector<vector<vector<string>>> &result, string exceptedCodons) {
        addBasesBase(std::move(aminoAcid), startElement, maxLength, result, BASE_EXTENDED, std::move(exceptedCodons));
    }

    void addOriginalBases(string aminoAcid, int startElement, unsigned int maxLength,
                          vector<vector<vector<string>>> &result, string exceptedCodons) {
        addBasesBase(std::move(aminoAcid), startElement, maxLength, result, BASE_ORIGINAL, std::move(exceptedCodons));
    }

    auto calcEachBaseEnergyBase(const string &aminoAcid, vector<vector<vector<string>>> &bases, int flag)
        -> vector<vector<vector<vector<pair<int, string>>>>> {
        // 結果格納ベクターを初期化
        int baseLength = aminoAcid.size() * 3;
        int maxLength = 8;
        vector<vector<vector<vector<pair<int, string>>>>> result;
        initResultExtendedBaseEnergyVector(baseLength, maxLength, result);

        // 各ポジション（配列長=5, 6, 8）における部分配列の自由エネルギーを取得
        vector<int> baseLengths;
        baseLengths.push_back(5);
        baseLengths.push_back(6);
        baseLengths.push_back(8);
        for (int baseSize : baseLengths) {
            for (int i = 0; i < baseLength; i++) {
                int position = i + 1;
                for (unsigned int k = 0; k < bases[position][baseSize].size(); k++) {
                    string base = bases[position][baseSize][k];

                    // 配列、エネルギーマップから自由エネルギーを検索
                    string originalBase = base;
                    if (flag == BASE_EXTENDED) {
                        // 拡張した塩基を使用する場合は、V、WはUに変換、X、YはGに変換
                        Util::baseReplace(originalBase, "V", "U");
                        Util::baseReplace(originalBase, "W", "U");
                        Util::baseReplace(originalBase, "X", "G");
                        Util::baseReplace(originalBase, "Y", "G");
                    }

                    // エネルギー取得（未設定なら以降の処理スキップ）
                    map<string, int>::iterator itr;
                    itr = baseEnergy.find(originalBase);
                    if (itr == baseEnergy.end()) {
                        continue;
                    }
                    int energy = itr->second;

                    // 結果ベクターに格納
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

    auto calcEachExtendedBaseEnergy(string aminoAcid, vector<vector<vector<string>>> &bases)
        -> vector<vector<vector<vector<pair<int, string>>>>> {
        return calcEachBaseEnergyBase(std::move(aminoAcid), bases, BASE_EXTENDED);
    }

    auto calcEachOriginalBaseEnergy(string aminoAcid, vector<vector<vector<string>>> &bases)
        -> vector<vector<vector<vector<pair<int, string>>>>> {
        return calcEachBaseEnergyBase(std::move(aminoAcid), bases, BASE_ORIGINAL);
    }

    auto getSelectedLengthExtendedBases(string aminoAcid, unsigned int maxLength, string exceptedCodons)
        -> vector<vector<vector<string>>> {
        vector<vector<vector<string>>> result;
        getBasesBase(std::move(aminoAcid), maxLength, result, BASE_EXTENDED, std::move(exceptedCodons));
        return result;
    }

    auto getSelectedLengthOriginalBases(string aminoAcid, unsigned int maxLength, string exceptedCodons)
        -> vector<vector<vector<string>>> {
        vector<vector<vector<string>>> result;
        getBasesBase(std::move(aminoAcid), maxLength, result, BASE_ORIGINAL, std::move(exceptedCodons));
        return result;
    }

    /*
     * デバッグ用
     */
    void printTime(const string &msg, clock_t start, clock_t end) {
        double margin = (double)(end - start) / CLOCKS_PER_SEC;
        cout << msg << "\t" << margin << endl;
    }
};

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

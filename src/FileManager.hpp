/*
 * FileManager.h
 *
 *  Created on: 2014/09/04
 *      Author: kamegai
 */

#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <regex.h>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>

using namespace std;

class FileManager {
  public:
    FileManager() = default;

    ~FileManager() = default;

    void saveSingleFastaFile(const string& name, const string& seq, const string& file) {
        ofstream ofs(file.c_str());
        ofs << ">" << name << endl;
        ofs << seq << endl;
    }

    auto loadFastaFile(const string& file) -> map<string, string>;

    void loadEnergyFile(map<string, int> &result);

  private:
    auto split(string &str, char delim) -> vector<string> {
        istringstream iss(str);
        string tmp;
        vector<string> res;

        while (getline(iss, tmp, delim))
            res.push_back(tmp);
        return res;
    }
};

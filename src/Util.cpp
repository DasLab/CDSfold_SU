/*
 * Util.cpp
 *
 *  Created on: 2021/2/22
 *      Author: everyday847
 */

#include <iostream>
#include <map>
#include <string>

#include "Util.hpp"

using namespace std;

void Util::baseReplace(string &base, const string& from, const string& to) {
    unsigned long int pos = 0;
    while ((pos = base.find(from, pos)) != string::npos) {
        base.replace(pos, from.length(), to);
        pos += to.length();
    }
}

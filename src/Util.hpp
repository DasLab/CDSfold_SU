/*
 * Util.h
 *
 *  Created on: 2014/09/22
 *      Author: kamegai
 */

#pragma once

#include <iostream>
#include <map>
#include <string>

using namespace std;

class Util {
  public:
    static void baseReplace(string &base, const string& from, const string& to) {
        int pos = 0;
        while ((pos = base.find(from, pos)) != string::npos) {
            base.replace(pos, from.length(), to);
            pos += to.length();
        }
    }
};

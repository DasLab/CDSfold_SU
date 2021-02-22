/*
 * Util.h
 *
 *  Created on: 2014/09/22
 *      Author: kamegai
 */

#pragma once

#include <iosfwd>


class Util {
    using string = std::string;
    public:
        static void baseReplace(string &base, const string& from, const string& to);
};

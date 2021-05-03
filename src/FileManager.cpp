#include "FileManager.hpp"


/* return the energy associated with various sequences */
void FileManager::loadEnergyFile(map<string, int> &result) {
    // Source: http://hackage.haskell.org/package/BiobaseViennaenergy_par.c
    
    /* Triloop37 */
    /* Number of loop partial bases (including base pair pairs), Base sequence free energy */

    result.insert(make_pair("CAACG", 680));
    result.insert(make_pair("GUUAC", 690));

    /* 6loop37 */
    result.insert(make_pair("CAACGG", 550));
    result.insert(make_pair("CCAAGG", 330));
    result.insert(make_pair("CCACGG", 370));
    result.insert(make_pair("CCCAGG", 340));
    result.insert(make_pair("CCGAGG", 350));
    result.insert(make_pair("CCGCGG", 360));
    result.insert(make_pair("CCUAGG", 370));
    result.insert(make_pair("CCUCGG", 250));
    result.insert(make_pair("CUAAGG", 360));
    result.insert(make_pair("CUACGG", 280));
    result.insert(make_pair("CUCAGG", 370));
    result.insert(make_pair("CUCCGG", 270));
    result.insert(make_pair("CUGCGG", 280));
    result.insert(make_pair("CUUAGG", 350));
    result.insert(make_pair("CUUCGG", 370));
    result.insert(make_pair("CUUUGG", 370));
    
    /* Hexaloop37 */
    result.insert(make_pair("ACAGUACU", 280));
    result.insert(make_pair("ACAGUGAU", 360));
    result.insert(make_pair("ACAGUGCU", 290));
    result.insert(make_pair("ACAGUGUU", 180));

}

/* load a fasta file and store it as a map from names to bases */
auto FileManager::loadFastaFile(const string& file) -> map<string, string> {
    string name = "";
    string base = "";
    map<string, string> fastaMap;
    ifstream ifs(file.c_str());

    /* get name from file stream and place into name */ 
    getline(ifs, name);

    /* regex pattern matching */ 
    regex_t preg;
    if (regcomp(&preg, "^>", REG_EXTENDED | REG_NEWLINE) != 0) {
        printf("regex compile failed.\n");
        // TODO Error handling
    }

    size_t nmatch = 1; // Maximum number of matches is one
    std::vector<regmatch_t> pmatch;
		pmatch.reserve(nmatch);
    
    string str;
    while (getline(ifs, str)) {
        if (regexec(&preg, str.c_str(), nmatch, &pmatch[0], 0) == 0) {
            fastaMap[name] = base;
            name = str;
            base = "";
        } else {
            base += str;
        }
    }

    fastaMap[name] = base;

    regfree(&preg);
    return fastaMap;
}

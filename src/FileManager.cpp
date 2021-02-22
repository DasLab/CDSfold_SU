#include "FileManager.hpp"


void FileManager::loadEnergyFile(map<string, int> &result) {
    /////////////////////////////////////////////////////////////////
    // 出典元 : http://hackage.haskell.org/package/BiobaseViennaenergy_par.c
    /////////////////////////////////////////////////////////////////
    //# Triloop37
    // ループ部分塩基数	ループ部分塩基数（塩基対ペア含む）	塩基配列	自由エネルギー
    result.insert(make_pair("CAACG", 680));
    result.insert(make_pair("GUUAC", 690));

    // 6loop37
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
    // Hexaloop37
    result.insert(make_pair("ACAGUACU", 280));
    result.insert(make_pair("ACAGUGAU", 360));
    result.insert(make_pair("ACAGUGCU", 290));
    result.insert(make_pair("ACAGUGUU", 180));
    //		// テスト用
    //		result.insert(make_pair("AUGCU", 505));
    //		result.insert(make_pair("AUGUU", 333));
    //		result.insert(make_pair("AUUUU", 111));
    //		result.insert(make_pair("UAUAUUUU", 222));
}

auto FileManager::loadFastaFile(const string& file) -> map<string, string> {
    string name = "";
    string base = "";
    map<string, string> fastaMap;
    ifstream ifs(file.c_str());

    // 最初のコンティグ名を取得
    getline(ifs, name);

    // パターンマッチ
    regex_t preg;
    if (regcomp(&preg, "^>", REG_EXTENDED | REG_NEWLINE) != 0) {
        printf("regex compile failed.\n");
        // TODO エラー処理
    }

    size_t nmatch = 1; // 「先頭>」のマッチ数は最大1
    std::vector<regmatch_t> pmatch;
		pmatch.reserve(nmatch);
    // unsigned int i; //, j;
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

    // 最後のコンティグをハッシュに設定
    fastaMap[name] = base;

    regfree(&preg);
    return fastaMap;
}
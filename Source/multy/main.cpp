#include "Aligners.h"
#include <iostream>

int main(int argc, char** argv)
{
    std::string str1 = "MGLLIALALLCLFSLAEANSKAITTSLTTKWFSAPLLLEASEFLAEDSQEKFWSFVEASQ";
    std::string str2 = "MLRAVALCVSVVLIALYTPTSGESSQSYPITTLINAKWTQTPLYLEIAEYLADEQAGLFW";
    NeedlmanWunsch nw(str1, str2, "BLOSUM62", -2, -1);
    string_tuple result = nw.Align();
    std::cout << result.first << std::endl;
    std::cout << result.second << std::endl;
    std::cout << nw.GetScore() << std::endl;
    return 0;
}

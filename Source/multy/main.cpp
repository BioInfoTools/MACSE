#include "Aligners.h"
#include <iostream>

int main(int argc, char** argv)
{
    std::string str1 = "ATCGAGATG";
    std::string str2 = "ATTTCGAAATG";
    NeedlmanWunsch nw(str1, str2, "BLOSUM62", -2, -1);
    Hein h(str1, str2, "NUC-45", "BLOSUM62", -2, -1, -2);
    //тестирование Нидлмана-Вунша
    string_tuple result = nw.Align();
    std::cout << "Needlman-Wunsch align:" << std::endl;
    std::cout << result.first << std::endl;
    std::cout << result.second << std::endl;
    //тестирование алгоритма Хейна
    result = h.Align();
    std::cout << "Hein align:" << std::endl;
    std::cout << result.first << std::endl;
    std::cout << result.second << std::endl;
    return 0;
}

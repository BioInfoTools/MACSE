#include "Aligners.h"

#include <iostream>
#include <stdlib.h>

int main(int argc, char** argv)
{
    if (argc == 9) {
        std::string str1(argv[1]);//"ATCGAGATG";
        std::string str2(argv[2]);//"ATTTCGAAATG";
        std::string NTsubs(argv[3]);
        std::string AAsubs(argv[4]);
        int gap_open = atoi(argv[5]);
        int gap_extension = atoi(argv[6]);
        int gap_frame = atoi(argv[7]);
        int stop_cost = atoi(argv[8]);
        //============================
        NeedlmanWunsch nw(str1, str2, NTsubs, gap_open, gap_extension);
        Hein h(str1, str2, NTsubs, AAsubs, stop_cost, gap_open, gap_extension, gap_frame);
        // тестирование Нидлмана-Вунша
        string_tuple result = nw.Align();
        std::cout << "Needlman-Wunsch align:" << std::endl;
        std::cout << result.first << std::endl;
        std::cout << result.second << std::endl;
        //тестирование алгоритма Хейна
        result = h.Align_old();
        std::cout << "Hein align (my align):" << std::endl;
        std::cout << result.first << std::endl;
        std::cout << result.second << std::endl;
        result = h.GetAAalign();
        std::cout << result.first << std::endl;
        std::cout << result.second << std::endl;
        //============================
        h.ChangeStrings(str1, str2);
        result = h.Align();
        std::cout << "Hein align (MACSE):" << std::endl;
        std::cout << result.first << std::endl;
        std::cout << result.second << std::endl;
        result = h.GetAAalign();
        std::cout << result.first << std::endl;
        std::cout << result.second << std::endl;
    } else {
        std::cout << "Wrong argument count!" << std::endl;
        std::cout << "./multy seq1 seq2 NTsubs AAsubs gap_open gap_extension gap_frame stop_cost" << std::endl;
    }
    return 0;
}

#include "Aligners.h"

#include <iostream>
#include <stdlib.h>

int main(int argc, char** argv) {
	if (argc > 2) {
		//имя программы + 2 строки на выравнивание
		std::string str1(argv[1]);
		std::string str2(argv[2]);
		std::string NTsubs = "NUC-45", AAsubs = "BLOSUM62";
		int gap_open = -2, gap_extension = -1, gap_frame = -2, stop_cost = -3;
		if (argc == 9) {
			//две матрицы замен + четыре штрафа
			NTsubs = argv[3];
			AAsubs = argv[4];
			gap_open = atoi(argv[5]);
			gap_extension = atoi(argv[6]);
			gap_frame = atoi(argv[7]);
			stop_cost = atoi(argv[8]);
		}
		//Нидлман-Вунш
		NeedlmanWunsch nw(str1, str2, NTsubs, gap_open, gap_extension);
		string_tuple result = nw.Align();
		std::cout << result.first  << std::endl;
		std::cout << result.second << std::endl;
		std::cout << nw.GetScore() << std::endl;
	} else std::cout << "Wrong argument count!" << std::endl;
	return 0;
}
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
		//Выравнивание
		string_tuple result;
		//Нидлман-Вунш
		/*
		NeedlmanWunsch nw(str1, str2, NTsubs, gap_open, gap_extension);
		nw.Align();
		result = nw.GetAAalign(AAsubs, gap_frame);
		*/
		Hein h(str1, str2, NTsubs, AAsubs, stop_cost, gap_open, gap_extension, gap_frame);
		//Hein 
		/*
		h.LoadHein();
		h.Align_old();
		result = h.GetAAalign();
		*/
		//MACSE
		h.LoadMACSE();
		h.Align();
		result = h.GetAAalign();
		//Ответ
		std::cout << result.first  << std::endl;
		std::cout << result.second << std::endl;
		std::cout << h.GetScore() << std::endl;
	} else std::cout << "Wrong argument count!" << std::endl;
	return 0;
}
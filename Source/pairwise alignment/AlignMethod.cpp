#include "Aligners.h"

#include <fstream>
#include <iostream>
#include <string.h>
#include <algorithm>

#define MAXINT 0x7FFFFFFF

void AlignMethod :: NewScoreMatrix(std::string file_name, int* score_matrix) {
	std::ifstream ifs;
	ifs.open (file_name, std::ifstream::in);
	
	//инициализация
	for (int i = 0; i < 128; i++)
		for (int j = 0; j < 128; j++)
			score_matrix[i*128 + j] = -MAXINT;
		
		//пропуск коментариев
		char c;
	std::string comment_string;
	while (ifs.good() && (c = ifs.get()) == '#')
		std::getline(ifs, comment_string);
	ifs.putback(c);
	
	//строка алфавита
	std::string alphabet;
	std::getline(ifs, alphabet);
	alphabet.erase(std::remove_if(alphabet.begin(), alphabet.end(),
								[](char c){ return (c == ' ' || c == '\t'); }), alphabet.end());
	
	//чтение таблицы
	for (unsigned int i = 0; i < alphabet.length(); i++) {
		c = ifs.get();
		if (c != '*')
			for (unsigned int j = 0; j < alphabet.length(); j++) {
				if (alphabet[j] != '*') ifs >> score_matrix[c*128 + alphabet[j]];
				else {
					int value;
					ifs >> value;
					for (int z = 0; z < 128; z++)
						if (score_matrix[c*128 + z] == -MAXINT)
							score_matrix[c*128 + z] = value;
				}
			}
			else
				for (unsigned int j = 0; j < alphabet.length(); j++) {
					int value;
					ifs >> value;
					if (alphabet[j] != '*') {
						for (int z = 0; z < 128; z++)
							if (score_matrix[z*128 + alphabet[j]] == -MAXINT)
								score_matrix[z*128 + alphabet[j]] = value;
					}
					else {
						for (int z = 0; z < 128; z++)
							if (score_matrix[z*128 + z] == -MAXINT)
								score_matrix[z*128 + z] = value;
					}
				}
				std::getline(ifs, comment_string);
	}
	
	ifs.close();
	/*  DEBUG
	 *    for (int i = 0; i < 128; i++) {
	 *        for (int j = 0; j < 128; j++) {
	 *            std::cout << score_matrix[i*128 + j] << ' ';
}
std::cout << std::endl;
}
*/
}

void AlignMethod :: IniFW() {
	memset(W, 0, sizeof(int)*n*m);
	memset(F, 0, sizeof(int)*n*m);
	/*
	 *    F[0] = gap_open;
	 *	for (int j = 1; j < m; j++) {
	 *        F[j] = F[j - 1] + gap_extension;
	 *        W[j] = 3;
}
for (int i = 1; i < n; i++) {
	F[i * m] = F[(i-1) * m] + gap_extension;
	W[i * m] = 2;
}
F[0] = 0;
*/
}

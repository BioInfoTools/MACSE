#include "Aligners.h"
#include <algorithm>
#include <string.h>
#include <fstream>
#include <iostream>


#define MAXINT 0x7FFFFFFF

void AlignMethod :: NewScoreMatrix(std::string file_name) {
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

    for (int i = 0; i < 128; i++) {
        for (int j = 0; j < 128; j++) {
            std::cout << score_matrix[i*128 + j] << ' ';
        }
        std::cout << std::endl;
    }
}

NeedlmanWunsch :: NeedlmanWunsch(std::string s1, std::string s2,
                                 std::string score_matrix,
                                 int gap_open, int gap_extension) {
    this->gap_open = gap_open;
	this->gap_extension = gap_extension;
	NewStrings(s1, s2);
    NewScoreMatrix(score_matrix);
}

string_tuple NeedlmanWunsch :: Align() {
    //прямой проход
    for (int i = 1; i < n; i++) {
		for (int j = 1; j < m; j++) {
			//оставить все как есть
			int way1 = F[(i-1)*m + j-1] +
                score_matrix[seq1[i-1]*128 + seq2[j-1]];
			//порвать seq1
			int way2 = F[(i-1)*m + j] + gap_extension;
			if (W[(i-1)*m + j] == 1) way2 += gap_open;
			//порвать seq2
			int way3 = F[i*m + j-1] + gap_extension;
			if (W[i*m + j-1] == 1) way3 += gap_open;
			//выбираем максимум
			if (way1 >= way2 && way1 >= way3) {
				F[i*m + j] = way1;
				W[i*m + j] = 1;
			} else if (way2 >= way1 && way2 >= way3) {
				F[i*m + j] = way2;
				W[i*m + j] = 2;
			} else {
				F[i*m + j] = way3;
				W[i*m + j] = 3;
			}
		}
	}
	//обратный ход, получение ответа
	result_score = F[(n-1)*m + m-1];
	int i = n-1, j = m-1;
	while (i > 0 || j > 0) {
        switch (W[i*m + j]) {
            case 1:
                align1 += seq1[--i];
                align2 += seq2[--j];
                break;
            case 2:
                align1 += seq1[--i];
                align2 += '-';
                break;
            case 3:
                align1 += '-';
                align2 += seq2[--j];
                break;
        }
    }
    std::reverse(align1.begin(), align1.end());
    std::reverse(align2.begin(), align2.end());
    return std::make_pair(align1, align2);
}

void NeedlmanWunsch :: NewStrings(std::string s1, std::string s2) {
    seq1 = s1; seq2 = s2;
    align1 = ""; align2 = "";
    result_score = 0;
    n = s1.length() + 1;
    m = s2.length() + 1;
    if (F) delete F;
    if (W) delete W;
    F = new int [n * m];
    W = new int [n * m];
    memset(W, 0, sizeof(int)*n*m);
    F[0] = gap_open;
	for (int j = 1; j < m; j++) {
        F[j] = F[j - 1] + gap_extension;
        W[j] = 3;
    }
	for (int i = 1; i < n; i++) {
        F[i * m] = F[(i-1) * m] + gap_extension;
        W[i * m] = 3;
    }
    F[0] = 0;
}

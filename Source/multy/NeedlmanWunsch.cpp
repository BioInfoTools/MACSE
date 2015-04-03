#include "Aligners.h"

#include <string.h>
#include <algorithm>

NeedlmanWunsch :: NeedlmanWunsch(std::string s1, std::string s2,
                                 std::string score_matrix,
                                 int gap_open, int gap_extension) {
    this->gap_open = gap_open;
	this->gap_extension = gap_extension;
	ChangeStrings(s1, s2);
    NewScoreMatrix(score_matrix, this->score_matrix);
}

string_tuple NeedlmanWunsch :: Align() {
    //прямой проход
    for (int i = 1; i < n; i++) {
		for (int j = 1; j < m; j++) {
			//оставить все как есть
			int way1 = F[(i-1)*m + j-1] +
                score_matrix[seq1[i-1]*128 + seq2[j-1]];
			//порвать seq2
			int way2 = F[(i-1)*m + j] + gap_extension;
			if (W[(i-1)*m + j] == 1) way2 += gap_open;
			//порвать seq1
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
	result_score = F[n*m-1];
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

void NeedlmanWunsch :: ChangeStrings(std::string s1, std::string s2) {
    seq1 = s1; seq2 = s2;
    align1 = ""; align2 = "";
    result_score = 0;
    n = s1.length() + 1;
    m = s2.length() + 1;
    if (F) delete F;
    if (W) delete W;
    F = new int [n * m];
    W = new int [n * m];
    IniFW();
}

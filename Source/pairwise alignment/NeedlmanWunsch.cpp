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

char NeedlmanWunsch :: TranslateNTtoAA(std::string& s, int i) {
		if (s[i] == '-' && s[i+1] == '-' && s[i+2] == '-') return '-';
		if (s[i] == '-' || s[i+1] == '-' || s[i+2] == '-') return '!';
    char result;
					switch (s[i]) {
            case 'A':
                switch (s[i+1]) {
                    case 'A':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result = 'K';
                                break;
                            default:
                                result = 'N';
                        }
                        break;
                    case 'C':
                        result = 'T';
                        break;
                    case 'G':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result = 'R';
                                break;
                            default:
                                result = 'S';
                        }
                        break;
                    default:
                        switch (s[i+2]) {
                            case 'G':
                                result = 'M'; //start
                                break;
                            default:
                                result = 'I';
                        }
                }
                break;
            case 'C':
                switch (s[i+1]) {
                    case 'A':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result = 'Q';
                                break;
                            default:
                                result = 'H';
                        }
                        break;
                    case 'C':
                        result = 'P';
                        break;
                    case 'G':
                        result = 'R';
                        break;
                    default:
                        result = 'L';
                }
                break;
            case 'G':
                switch (s[i+1]) {
                    case 'A':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result = 'E';
                                break;
                            default:
                                result = 'D';
                        }
                        break;
                    case 'C':
                        result = 'A';
                        break;
                    case 'G':
                        result = 'G';
                        break;
                    default:
                        result = 'V';
                }
                break;
            default:
                switch (s[i+1]) {
                    case 'A':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result = '*'; //stop
                                break;
                            default:
                                result = 'Y';
                        }
                        break;
                    case 'C':
                        result = 'S';
                        break;
                    case 'G':
                        switch (s[i+2]) {
                            case 'A':
                                result = '*'; //stop
                                break;
                            case 'G':
                                result = 'W';
                                break;
                            default:
                                result = 'C';
                        }
                        break;
                    default:
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result = 'L';
                                break;
                            default:
                                result = 'F';
                        }
                }
        }
   
    return result;
}


string_tuple NeedlmanWunsch :: GetAAalign(std::string AAscore_matrix, int gap) {
	//загружаем матрицу замен 
	int aa_score_matrix[128*128];
	NewScoreMatrix(AAscore_matrix, aa_score_matrix);
	//выбор рамки
	int max_score = 0x80000000, shift;
	for (int frame_shift = 0; frame_shift < 3; frame_shift++) {
			int tek_score = 0;
			for (int i = frame_shift; i < align1.length() - 2; i += 3) {
					char aa1 = TranslateNTtoAA(align1, i);
					char aa2 = TranslateNTtoAA(align2, i);
					if (aa1 == '!') tek_score += gap;
					if (aa2 == '!') tek_score += gap;
					if (aa1 == '-') tek_score += gap_extension;
					if (aa2 == '-') tek_score += gap_extension;
					if (aa1 != '!' && aa1 != '-' && aa2 != '!' && aa2 != '-')
							tek_score += aa_score_matrix[aa1*128 + aa2];
			}
			if (tek_score > max_score) {
				max_score = tek_score;
				shift = frame_shift;
			}
	}
	std::string result1 = "";
	std::string result2 = "";
	if (shift) {
		result1 += '!';
		result2 += '!';
	}
	for (int i = shift; i < align1.length() - 2; i += 3) {
			result1 += TranslateNTtoAA(align1, i);
			result2 += TranslateNTtoAA(align2, i);
	}
	if ((align1.length() - shift) % 3) {
		result1 += '!';
		result2 += '!';
	}
	return std::make_pair(result1, result2);
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

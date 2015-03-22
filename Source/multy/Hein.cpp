#include "Aligners.h"

#include <string.h>
#include <algorithm>

Hein :: Hein(std::string s1, std::string s2,
             std::string nt_score_matrix, std::string aa_score_matrix,
             int gap_open, int gap_extension, int gap_frame) {
    this->gap_open = gap_open;
    this->gap_frame = gap_frame;
	this->gap_extension = gap_extension;
	ChangeStrings(s1, s2);
    ChangeNTscoreMatrix(nt_score_matrix);
    ChangeAAscoreMatrix(aa_score_matrix);
}

string_tuple Hein :: Align() {
    for (int i = 1; i < 3 && i < n; i++)
        for (int j = 1; j < m; j++) {
            //оставить все как есть
			int way5 = F[(i-1)*m + j-1] + gap_frame +
                nt_score_matrix[seq1[i-1]*128 + seq2[j-1]];
			//порвать seq2
			int way4 = F[(i-1)*m + j] + gap_extension + gap_frame;
			if (W[(i-1)*m + j] > 4) way4 += gap_open;
			//порвать seq1
			int way3 = F[i*m + j-1] + gap_extension + gap_frame;
			if (W[i*m + j-1] > 4) way3 += gap_open;
			//выбираем максимум
			if (way5 >= way4 && way5 >= way3) {
				F[i*m + j] = way5;
				W[i*m + j] = 5;
			} else if (way4 >= way5 && way4 >= way3) {
				F[i*m + j] = way4;
				W[i*m + j] = 4;
			} else {
				F[i*m + j] = way3;
				W[i*m + j] = 3;
			}
        }
    for (int i = 3; i < n; i++)
        for (int j = 1; j < 3 && j < m; j++) {
            //оставить все как есть
			int way5 = F[(i-1)*m + j-1] + gap_frame +
                nt_score_matrix[seq1[i-1]*128 + seq2[j-1]];
			//порвать seq2
			int way4 = F[(i-1)*m + j] + gap_extension + gap_frame;
			if (W[(i-1)*m + j] > 4) way4 += gap_open;
			//порвать seq1
			int way3 = F[i*m + j-1] + gap_extension + gap_frame;
			if (W[i*m + j-1] > 4) way3 += gap_open;
			//выбираем максимум
			if (way5 >= way4 && way5 >= way3) {
				F[i*m + j] = way5;
				W[i*m + j] = 5;
			} else if (way4 >= way5 && way4 >= way3) {
				F[i*m + j] = way4;
				W[i*m + j] = 4;
			} else {
				F[i*m + j] = way3;
				W[i*m + j] = 3;
			}
        }
    //прямой проход
    for (int i = 3; i < n; i++) {
		for (int j = 3; j < m; j++) {
			//AA gap seq1
            int max = F[i*m + j-3] + gap_extension*3, index = 1;
            if (W[i*m + j-3] > 4) max += gap_open;
            //AA gap seq2
            int score = F[(i-3)*m + j] + gap_extension*3;
            if (W[(i-3)*m + j] > 4) score += gap_open;
            if (score > max) {
                max = score;
                index = 2;
            }
            //NT gap seq1
            score = F[i*m + j-1] + gap_extension + gap_frame;
            if (W[i*m + j-1] > 4) score += gap_open;
            if (score > max) {
                max = score;
                index = 3;
            }
            //NT gap seq2
            score = F[(i-1)*m + j] + gap_extension + gap_frame;
            if (W[(i-1)*m + j] > 4) score += gap_open;
            if (score > max) {
                max = score;
                index = 4;
            }
            //NT совпадение
            score = F[(i-1)*m + j-1] + gap_frame +
                nt_score_matrix[seq1[i-1]*128 + seq2[j-1]];
            if (W[(i-1)*m + j-1] > 4) score += gap_open;
            if (score > max) {
                max = score;
                index = 5;
            }
            //AA совпадение
            score = F[(i-3)*m + j-3] +
                aa_score_matrix[AAseq1[i-3]*128 + AAseq2[j-3]] +
                nt_score_matrix[seq1[i-3]*128 + seq2[j-3]] +
                nt_score_matrix[seq1[i-2]*128 + seq2[j-2]] +
                nt_score_matrix[seq1[i-1]*128 + seq2[j-1]];
            if (W[(i-1)*m + j-1] > 4) score += gap_open;
            if (score > max) {
                max = score;
                index = 6;
            }
            //сохранение результата
            F[i*m + j] = max;
            W[i*m + j] = index;
		}
	}
	//обратный ход, получение ответа
	result_score = F[n*m-1];
	int i = n-1, j = m-1;
	while (i > 0 || j > 0) {
        switch (W[i*m + j]) {
            case 1:
                //AA gap seq1
                align1 += "---";
                align2 += seq2[--j];
                align2 += seq2[--j];
                align2 += seq2[--j];
                break;
            case 2:
                //AA gap seq2
                align1 += seq1[--i];
                align1 += seq1[--i];
                align1 += seq1[--i];
                align2 += "---";
                break;
            case 3:
                //NT gap seq1
                align1 += '-';
                align2 += seq1[--j];
                break;
            case 4:
                //NT gap seq2
                align1 += seq1[--i];
                align2 += '-';
                break;
            case 5:
                //NT match
                align1 += seq1[--i];
                align2 += seq2[--j];
                break;
            case 6:
                //AA match
                align1 += seq1[--i];
                align1 += seq1[--i];
                align1 += seq1[--i];
                align2 += seq2[--j];
                align2 += seq2[--j];
                align2 += seq2[--j];
                break;
        }
    }
    std::reverse(align1.begin(), align1.end());
    std::reverse(align2.begin(), align2.end());

    return std::make_pair(align1, align2);
}

//пздц... <(_ _)>
std::string Hein :: TranslateNTtoAA(std::string& s) {
    std::string result = "";
    for (unsigned int i = 0; i < s.length() - 2; i++) {
        switch (s[i]) {
            case 'A':
                switch (s[i+1]) {
                    case 'A':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result += 'K';
                                break;
                            default:
                                result += 'N';
                        }
                        break;
                    case 'C':
                        result += 'T';
                        break;
                    case 'G':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result += 'R';
                                break;
                            default:
                                result += 'S';
                        }
                        break;
                    default:
                        switch (s[i+2]) {
                            case 'G':
                                result += 'M'; //start
                                break;
                            default:
                                result += 'I';
                        }
                }
                break;
            case 'C':
                switch (s[i+1]) {
                    case 'A':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result += 'Q';
                                break;
                            default:
                                result += 'H';
                        }
                        break;
                    case 'C':
                        result += 'P';
                        break;
                    case 'G':
                        result += 'R';
                        break;
                    default:
                        result += 'L';
                }
                break;
            case 'G':
                switch (s[i+1]) {
                    case 'A':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result += 'E';
                                break;
                            default:
                                result += 'D';
                        }
                        break;
                    case 'C':
                        result += 'A';
                        break;
                    case 'G':
                        result += 'G';
                        break;
                    default:
                        result += 'V';
                }
                break;
            default:
                switch (s[i+1]) {
                    case 'A':
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result += '*'; //stop
                                break;
                            default:
                                result += 'Y';
                        }
                        break;
                    case 'C':
                        result += 'S';
                        break;
                    case 'G':
                        switch (s[i+2]) {
                            case 'A':
                                result += '*'; //stop
                                break;
                            case 'G':
                                result += 'W';
                                break;
                            default:
                                result += 'C';
                        }
                        break;
                    default:
                        switch (s[i+2]) {
                            case 'A':
                            case 'G':
                                result += 'L';
                                break;
                            default:
                                result += 'F';
                        }
                }
        }
    }
    return result;
}

void Hein :: ChangeStrings(std::string s1, std::string s2) {
    seq1 = s1; seq2 = s2;
    align1 = ""; align2 = "";
    AAseq1 = TranslateNTtoAA(seq1);
    AAseq2 = TranslateNTtoAA(seq2);
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
        W[i * m] = 2;
    }
    F[0] = 0;
    W[0] = 0;
}

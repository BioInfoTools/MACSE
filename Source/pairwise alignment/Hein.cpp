#include "Aligners.h"

#include <string.h>
#include <algorithm>

#define MAXINT 0x7FFFFFFF

Hein :: Hein(std::string s1, std::string s2,
	std::string nt_score_matrix, std::string aa_score_matrix,
	int stop_cost, int gap_open, int gap_extension, int gap_frame) {
	this->gap_open = gap_open;
	this->gap_frame = gap_frame;
	this->gap_extension = gap_extension;
	this->stop_cost = stop_cost;
	ChangeStrings(s1, s2);
	ChangeNTscoreMatrix(nt_score_matrix);
	ChangeAAscoreMatrix(aa_score_matrix);
}


//ебучий костыль #1
std::string Kostyl(std::string& s) {
	std::string result = "";
	for (unsigned int i = 0; i < s.length(); i += 3) {
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
	if (s.length() % 3) result += '!';
	return result;
}

//ебучий костыль #2
std::string Kostyl2(std::string s) {
	std::string result = "";
	for (int i = s.length()-1; i > 1; i -= 3) {
		switch (s[i]) {
			case 'A':
				switch (s[i-1]) {
					case 'A':
						switch (s[i-2]) {
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
								switch (s[i-2]) {
									case 'A':
									case 'G':
										result += 'R';
										break;
									default:
										result += 'S';
								}
								break;
									default:
										switch (s[i-2]) {
											case 'G':
												result += 'M'; //start
												break;
											default:
												result += 'I';
										}
				}
				break;
											case 'C':
												switch (s[i-1]) {
													case 'A':
														switch (s[i-2]) {
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
																switch (s[i-1]) {
																	case 'A':
																		switch (s[i-2]) {
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
																				switch (s[i-1]) {
																					case 'A':
																						switch (s[i-2]) {
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
																								switch (s[i-2]) {
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
																										switch (s[i-2]) {
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
	if (s.length() % 3) result += '!';
	std::reverse(result.begin(), result.end());
	return result;
}



string_tuple Hein :: Align() {
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < m; j++) {
			//трансляция
			char AA1 = '?', AA2 = '?';
			if (i > 2) AA1 = AAseq1[i-3];
			if (j > 2) AA2 = AAseq2[j-3];
			//проверка на стоп-кодоны
			int stopS1, stopS2;
			stopS1 = (AA1 == '*') ? stop_cost : 0;
			stopS2 = (AA2 == '*') ? stop_cost : 0;
			//замена на AA уровне
			int subst_AA;
			if (AA1 == '*' || AA2 == '*') subst_AA = stopS1 + stopS2;
			else if (i - 3 >= 0 && j - 3 >= 0)
				subst_AA = aa_score_matrix[AA1*128 + AA2];
			else subst_AA = MAXINT;
			//вычисляем оптимальный ход
			//NT align
			int score = F[(i-1)*m + j - 1] + nt_score_matrix[seq1[i-1]*128+seq2[j-1]];
			int way = 12;
			if (j - 2 >= 0 && score < F[(i-1)*m + j-2]) {
				score = F[(i-1)*m + j - 2];
				way = 13;
			}
			if (i - 2 >= 0 && score < F[(i-2)*m + j-1]) {
				score = F[(i-2)*m + j-1];
				way = 14;
			}
			if (i - 2 >= 0 && j-2 >= 0 && score < F[(i-2)*m + j-2]) {
				score = F[(i-2)*m + j-2];
				way = 15;
			}
			score += 2 * gap_frame;
			//AA align
			int tmp;
			if (i-3 >= 0) {
				tmp = F[(i-3)*m + j] + stopS1 + gap_extension;
				if (W[(i-3)*m + j] > 6) tmp += gap_open;
				if (score <= tmp) {
					score = tmp;
					way = 1;
				}
			}
			if (j-3 >= 0) {
				tmp = F[i*m + j-3] + stopS2 + gap_extension;
				if (W[i*m + j-3] > 6) tmp += gap_open;
				if (score <= tmp) {
					score = tmp;
					way = 2;
				}
			}
			if (i-3 >= 0 && j-3 >= 0 && score <= F[(i-3)*m+j-3] + subst_AA) {
				score = F[(i-3)*m+j-3] + subst_AA;
				way = 7;
			}
			//AA & NT
			if (i > 2 && j > 1 && score < F[(i-3)*m+j-2]+stopS1+gap_frame) {
				score = F[(i-3)*m+j-2]+stopS1+gap_frame;
				way = 8;
			}
			if (i > 2 && score < F[(i-3)*m + j-1] + stopS1 + gap_frame) {
				score = F[(i-3)*m + j-1] + stopS1 + gap_frame;
				way = 9;
			}
			if (i > 1 && j > 2 && score < F[(i-2)*m+j-3]+stopS2+gap_frame) {
				score = score < F[(i-2)*m+j-3]+stopS2+gap_frame;
				way = 10;
			}
			if (j > 2 && score < F[(i-1)*m + j-3] + stopS2 + gap_frame) {
				score = F[(i-1)*m+j-3]+stopS2+gap_frame;
				way = 11;
			}
			//NT gap
			tmp = F[i*m+j-1] + gap_extension + gap_frame;
			if (W[i*m+j-1] > 6) tmp += gap_open;
			if (score < tmp) {
				score = tmp;
				way = 3;
			}
			if (j > 1) {
				tmp = F[i*m+j-2] + gap_extension + gap_frame;
				if (W[i*m+j-2] > 6) tmp += gap_open;
				if (score < tmp) {
					score = tmp;
					way = 4;
				}
			}
			tmp = F[(i-1)*m+j] + gap_extension + gap_frame;
			if (W[(i-1)*m+j] > 6) tmp += gap_open;
			if (score < tmp) {
				score = tmp;
				way = 5;
			}
			if (i > 1) {
				tmp = F[(i-2)*m+j] + gap_extension + gap_frame;
				if (W[(i-2)*m+j] > 6) tmp += gap_open;
				if (score < tmp) {
					score = tmp;
					way = 6;
				}
			}
			//сохраняем оптимальный шаг
			F[i*m+j] = score;
			W[i*m+j] = way;
		}
	}
	//обратный ход, получение ответа
	result_score = F[n*m-1];
	int i = n-1, j = m-1;
	//поиск минимума в последней строке
	for (int index = 0; index < m-1; index++) 
		if (F[(n-1)*m+index] > result_score) {
			j = index;
			result_score = F[(n-1)*m+index];
		}
	//поиск минимума в последнем столбце
	for (int index = 0; index < n-1; index++) 
		if (F[index*m + m-1] > result_score) {
			i = index; 
			j = m-1;
			result_score = F[index*m + m-1];
		}
	//i, j - точка минимума
	if (m-1-j) {
		align1 = std::string(m-1-j, '-');
		if ((m-1-j)%3) AAalign1 = '-';
		AAalign1 += std::string((m-1-j)/3, '-');
		align2 = seq2.substr(j, seq2.length() - j);
		AAalign2 = Kostyl(align2);
		std::reverse(align2.begin(), align2.end());
		std::reverse(AAalign2.begin(), AAalign2.end());
	} else {
		align1 = seq1.substr(i, seq1.length() - i);
		AAalign1 = Kostyl(align1);
		std::reverse(align1.begin(), align1.end());
		std::reverse(AAalign1.begin(), AAalign1.end());
		align2 = std::string(n-1-i, '-');
		if ((n-1-i)%3) AAalign2 = '-';
		AAalign2 += std::string((n-1-i)/3, '-');
	}
	int count1, count2, count3;
	while (W[i*m+j]) {
		switch (W[i*m + j]) {
			//AA level
			case 7:
				//AA match
				align1 += seq1[--i];
				align1 += seq1[--i];
				align1 += seq1[--i];
				align2 += seq2[--j];
				align2 += seq2[--j];
				align2 += seq2[--j];
				AAalign1 += AAseq1[i];
				AAalign2 += AAseq2[j];
				break;
			case 1:
				//AA gap seq2
				align1 += seq1[--i];
				align1 += seq1[--i];
				align1 += seq1[--i];
				align2 += "---";
				AAalign1 += AAseq1[i];
				AAalign2 += '-';
				break;
			case 2:
				//AA gap seq1
				align2 += seq2[--j];
				align2 += seq2[--j];
				align2 += seq2[--j];
				align1 += "---";
				AAalign1 += '-';
				AAalign2 += AAseq2[j];
				break;
				//AA & NT
			case 8:
				//NT gap seq2
				align1 += seq1[--i];
				align1 += seq1[--i];
				align1 += seq1[--i];
				AAalign1 += AAseq1[i];
				AAalign2 += '!';
				count1 = nt_score_matrix[seq1[i]*128 + seq2[j-2]] +
				nt_score_matrix[seq1[i+1]*128 + seq2[j-1]];
				count2 = nt_score_matrix[seq1[i]*128 + seq2[j-2]] +
				nt_score_matrix[seq1[i+2]*128 + seq2[j-1]];
				count3 = nt_score_matrix[seq1[i+1]*128 + seq2[j-2]] +
				nt_score_matrix[seq1[i+2]*128 + seq2[j-1]];
				if (count1 > count2 && count1 > count3) {
					align2 += '-';
					align2 += seq2[--j];
					align2 += seq2[--j];
				} else if (count2 > count1 && count2 > count3) {
					align2 += seq2[--j];
					align2 += '-';
					align2 += seq2[--j];
				} else {
					align2 += seq2[--j];
					align2 += seq2[--j];
					align2 += '-';
				}
				break;
			case 9:
				//NT gap seq2
				align1 += seq1[--i];
				align1 += seq1[--i];
				align1 += seq1[--i];
				AAalign1 += AAseq1[i];
				AAalign2 += '!';
				count1 = nt_score_matrix[seq1[i]*128 + seq2[j-1]];
				count2 = nt_score_matrix[seq1[i+1]*128 + seq2[j-1]];
				count3 = nt_score_matrix[seq1[i+2]*128 + seq2[j-1]];
				if (count1 > count2 && count1 > count3) {
					align2 += '-';
					align2 += '-';
					align2 += seq2[--j];
				} else if (count2 > count1 && count2 > count3) {
					align2 += '-';
					align2 += seq2[--j];
					align2 += '-';
				} else {
					align2 += seq2[--j];
					align2 += '-';
					align2 += '-';
				}
				break;
			case 10:
				//NT gap seq1
				align2 += seq2[--j];
				align2 += seq2[--j];
				align2 += seq2[--j];
				AAalign2 += AAseq2[j];
				AAalign1 += '!';
				count1 = nt_score_matrix[seq1[i-2]*128 + seq2[j]] +
				nt_score_matrix[seq1[i-1]*128 + seq2[j+1]];
				count2 = nt_score_matrix[seq1[i-2]*128 + seq2[j]] +
				nt_score_matrix[seq1[i-1]*128 + seq2[j+2]];
				count3 = nt_score_matrix[seq1[i-2]*128 + seq2[j+1]] +
				nt_score_matrix[seq1[i-1]*128 + seq2[j+2]];
				if (count1 > count2 && count1 > count3) {
					align1 += '-';
					align1 += seq1[--i];
					align1 += seq1[--i];
				} else if (count2 > count1 && count2 > count3) {
					align1 += seq1[--i];
					align1 += '-';
					align1 += seq1[--i];
				} else {
					align1 += seq1[--i];
					align1 += seq1[--i];
					align1 += '-';
				}
				break;
			case 11:
				//NT gap seq1
				align2 += seq2[--j];
				align2 += seq2[--j];
				align2 += seq2[--j];
				AAalign2 += AAseq2[j];
				AAalign1 += '!';
				count1 = nt_score_matrix[seq1[i-1]*128 + seq2[j]];
				count2 = nt_score_matrix[seq1[i-1]*128 + seq2[j+1]];
				count3 = nt_score_matrix[seq1[i-1]*128 + seq2[j+2]];
				if (count1 > count2 && count1 > count3) {
					align1 += '-';
					align1 += '-';
					align1 += seq1[--i];
				} else if (count2 > count1 && count2 > count3) {
					align1 += '-';
					align1 += seq1[--i];
					align1 += '-';
				} else {
					align1 += seq1[--i];
					align1 += '-';
					align1 += '-';
				}
				break;
				//NT gap
			case 3:
				align1 += '-';
				align2 += seq2[--j];
				AAalign1 += '-';
				AAalign2 += '!';
				break;
			case 4:
				align1 += '-';
				align1 += '-';
				align2 += seq2[--j];
				align2 += seq2[--j];
				AAalign1 += '-';
				AAalign2 += '!';
				break;
			case 5:
				align1 += seq1[--i];
				align2 += '-';
				AAalign1 += '!';
				AAalign2 += '-';
				break;
			case 6:
				align1 += seq1[--i];
				align1 += seq1[--i];
				align2 += '-';
				align2 += '-';
				AAalign1 += '!';
				AAalign2 += '-';
				break;
				//NT level
			case 12:
				align1 += seq1[--i];
				align2 += seq2[--j];
				AAalign1 += '!';
				AAalign2 += '!';
				break;
			case 13:
				align2 += seq2[--j];
				align2 += seq2[--j];
				count1 = nt_score_matrix[seq1[i-1]*128 + seq2[j]];
				count2 = nt_score_matrix[seq1[i-1]*128 + seq2[j+1]];
				if (count1 > count2) {
					align1 += '-';
					align1 += seq1[--i];
				} else {
					align1 += seq1[--i];
					align1 += '-';
				}
				AAalign1 += '!';
				AAalign2 += '!';
				break;
			case 14:
				align1 += seq1[--i];
				align1 += seq1[--i];
				count1 = nt_score_matrix[seq1[i]*128 + seq2[j-1]];
				count2 = nt_score_matrix[seq1[i+1]*128 + seq2[j-1]];
				if (count1 > count2) {
					align2 += '-';
					align2 += seq2[--j];
				} else {
					align2 += seq2[--j];
					align2 += '-';
				}
				AAalign1 += '!';
				AAalign2 += '!';
				break;
			case 15:
				align1 += seq1[--i];
				align1 += seq1[--i];
				align2 += seq2[--j];
				align2 += seq2[--j];
				AAalign1 += '!';
				AAalign2 += '!';
				break;
		}
	}
	align1 += std::string(j, '-');
	if (j%3) AAalign1 += '-';
	AAalign1 += std::string(j/3, '-');
	align2 += std::string(i, '-');
	if (i%3) AAalign2 += '-';
	AAalign2 += std::string(i/3, '-');
	
	std::reverse(align1.begin(), align1.end());
	std::reverse(align2.begin(), align2.end());
	
	align1 = seq1.substr(0, i) + align1;
	align2 = seq2.substr(0, j) + align2;
	
	std::reverse(AAalign1.begin(), AAalign1.end());
	std::reverse(AAalign2.begin(), AAalign2.end());
	
	AAalign1 = Kostyl2(seq1.substr(0, i)) + AAalign1;
	AAalign2 = Kostyl2(seq2.substr(0, j)) + AAalign2;
	
	return std::make_pair(align1, align2);
}

string_tuple Hein :: Align_old() {
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
	//поиск минимума в последней строке
	for (int index = 0; index < m-1; index++) 
		if (F[(n-1)*m+index] > result_score) {
			j = index;
			result_score = F[(n-1)*m+index];
		}
	//поиск минимума в последнем столбце
	for (int index = 0; index < n-1; index++) 
		if (F[index*m + m-1] > result_score) {
			i = index; 
			j = m-1;
			result_score = F[index*m + m-1];
		}
	//i, j - точка минимума
	if (m-1-j) {
		align1 = std::string(m-1-j, '-');
		if ((m-1-j)%3) AAalign1 = '-';
		AAalign1 += std::string((m-1-j)/3, '-');
		align2 = seq2.substr(j, seq2.length() - j);
		AAalign2 = Kostyl(align2);
		std::reverse(align2.begin(), align2.end());
		std::reverse(AAalign2.begin(), AAalign2.end());
	} else {
		align1 = seq1.substr(i, seq1.length() - i);
		AAalign1 = Kostyl(align1);
		std::reverse(align1.begin(), align1.end());
		std::reverse(AAalign1.begin(), AAalign1.end());
		align2 = std::string(n-1-i, '-');
		if ((n-1-i)%3) AAalign2 = '-';
		AAalign2 += std::string((n-1-i)/3, '-');
	}
	int len_mark = 0;
	while (W[i*m+j]) {
		switch (W[i*m + j]) {
			case 1:
				//AA gap seq1
				align1 += "---";
				align2 += seq2[--j];
				align2 += seq2[--j];
				align2 += seq2[--j];
				AAalign1 += '-';
				AAalign2 += AAseq2[j];
				len_mark = 0;
				break;
			case 2:
				//AA gap seq2
				align1 += seq1[--i];
				align1 += seq1[--i];
				align1 += seq1[--i];
				align2 += "---";
				AAalign1 += AAseq1[i];
				AAalign2 += '-';
				len_mark = 0;
				break;
			case 3:
				//NT gap seq1
				align1 += '-';
				align2 += seq2[--j];
				if (!len_mark) {
					AAalign1 += '-';
					AAalign2 += '!';
				}
				len_mark = (len_mark + 1) % 3;
				break;
			case 4:
				//NT gap seq2
				align1 += seq1[--i];
				align2 += '-';
				if (!len_mark) {
					AAalign1 += '!';
					AAalign2 += '-';
				}
				len_mark = (len_mark + 1) % 3;
				break;
			case 5:
				//NT match
				align1 += seq1[--i];
				align2 += seq2[--j];
				if (!len_mark) {
					AAalign1 += '!';
					AAalign2 += '!';
				}
				len_mark = (len_mark + 1) % 3;
				break;
			case 6:
				//AA match
				align1 += seq1[--i];
				align1 += seq1[--i];
				align1 += seq1[--i];
				align2 += seq2[--j];
				align2 += seq2[--j];
				align2 += seq2[--j];
				AAalign1 += AAseq1[i];
				AAalign2 += AAseq2[j];
				len_mark = 0;
				break;
		}
	}
	align1 += std::string(j, '-');
	if (j%3) AAalign1 += '-';
	AAalign1 += std::string(j/3, '-');
	align2 += std::string(i, '-');
	if (i%3) AAalign2 += '-';
	AAalign2 += std::string(i/3, '-');
	
	std::reverse(align1.begin(), align1.end());
	std::reverse(align2.begin(), align2.end());
	
	align1 = seq1.substr(0, i) + align1;
	align2 = seq2.substr(0, j) + align2;
	
	std::reverse(AAalign1.begin(), AAalign1.end());
	std::reverse(AAalign2.begin(), AAalign2.end());
	
	AAalign1 = Kostyl2(seq1.substr(0, i)) + AAalign1;
	AAalign2 = Kostyl2(seq2.substr(0, j)) + AAalign2;
	
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
	align1 = align2 = "";
	AAalign1 = AAalign2 = "";
	AAseq1 = TranslateNTtoAA(seq1);
	AAseq2 = TranslateNTtoAA(seq2);
	result_score = 0;
	n = s1.length() + 1;
	m = s2.length() + 1;
	if (F) delete F;
	if (W) delete W;
	F = new int [n * m];
	W = new int [n * m];
}

void Hein :: LoadMACSE() {
	align1 = align2 = AAalign1 = AAalign2 = "";
	memset(W, 0, sizeof(int)*n*m);
	memset(F, 0, sizeof(int)*n*m);
	/*
	*	F[0] = gap_open;
	*	for (int i = 3; i < n; i++) {
	*		F[i * m] = F[(i-3) * m] + gap_extension;
	*		if (AAseq1[i-3] == '*') F[i * m] += stop_cost;
	*		W[i * m] = 1;
}
for (int j = 3; j < m; j++) {
	F[j] = F[j - 3] + gap_extension;
	if (AAseq2[j-3] == '*') F[j] += stop_cost;
	W[j] = 2;
}
if (m) {
	F[1] = gap_open + gap_extension + gap_frame;
	W[1] = 3;
	if (m > 1) {
		F[2] = gap_open + gap_extension + gap_frame;
		W[2] = 4; 
}
}
if (n) {
	F[m] = gap_open + gap_extension + gap_frame;
	W[m] = 5;
	if (n > 1) {
		F[2*m] = gap_open + gap_extension + gap_frame;
		W[2*m] = 6;
}
}
F[0] = 0;
*/
}

void Hein :: LoadHein() {
	align1 = align2 = AAalign1 = AAalign2 = "";
	memset(W, 0, sizeof(int)*n*m);
	memset(F, 0, sizeof(int)*n*m);
	/*
	*	F[0] = gap_open;
	*	for (int j = 3; j < m; j++) {
	*		F[j] = F[j - 3] + gap_extension;
	*		W[j] = 1;
}
for (int i = 3; i < n; i++) {
	F[i * m] = F[(i-3) * m] + gap_extension;
	W[i * m] = 2;
}	
if (n) {
	F[1] = gap_open + gap_extension + gap_frame;
	W[1] = 3;
	if (n > 1) {
		F[2] = gap_open + gap_extension + gap_frame;
		W[2] = 3;
}
}
if (m) {
	F[m] = gap_open + gap_extension + gap_frame;
	W[m] = 4;
	if (m > 1) {
		F[2*m] = gap_open + gap_extension + gap_frame;
		W[2*m] = 4;
}
}
F[0] = 0;
*/
}
						 
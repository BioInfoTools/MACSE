#include <string.h>

#include "Profile.h"
#include "PairwiseAlign.h"

#define MAXINT 0x7FFFFFFF

#define DEBUG

#ifdef DEBUG
	#include <iostream>
	using std::cout;
	using std::endl;
#endif

// global tables for merge profiles
// memory allocation in Profile::operator +
// free in UPGMAfree
unsigned int* best_mvt = NULL;
float* scores = NULL;
int dim = 0;

// array of parameters for the aligner profiles 
int parameters[4]; 

void UPGMA(std::vector<BioSeq*>& sequences, PairwiseAlign& aligner) {
	int* matrix = new int[sequences.size() * sequences.size()];
	Profile* profiles = new Profile[sequences.size()];
	int max = -MAXINT;
	unsigned int max_i, max_j;
	for (unsigned int i = 0; i < sequences.size(); i++) {
		// at first each sequence is a single cluster
		profiles[i].sequences.push_back(sequences[i]);
		profiles[i].nt_score_matrix = aligner.GetNTscoreMatrix();
		profiles[i].aa_score_matrix = aligner.GetAAscoreMatrix();
		// calculating distances between other sequences
		for (unsigned int j = 0; j < i; j++) {
			matrix[i*sequences.size() + j] = aligner.Align(sequences[i],sequences[j]);
			// searching max score
			if (matrix[i*sequences.size() + j] > max) {
				max = matrix[i*sequences.size() + j];
				max_i = i;
				max_j = j;
			}
		}
	}
	
#ifdef DEBUG
	cout << "Building UPGMA matrix done!" << endl;
	for (unsigned int i = 0; i < sequences.size(); i++) {
		for (unsigned int j = 0; j < i; j++) {
			cout << matrix[i*sequences.size() + j] << '\t';
		}
		cout << '.' << endl;
	}
	cout << "start point:" << endl;
	cout << "max_i = " << max_i << endl;
	cout << "max_j = " << max_j << endl;
#endif

	parameters[0] = aligner.GetGapOpen();
	parameters[1] = aligner.GetGapExtension();
	parameters[2] = aligner.GetGapFrame();
	parameters[3] = aligner.GetStopCost();
	int alive = sequences.size();
	bool* corpse = new bool[sequences.size()];
	memset(corpse, false, sizeof(bool)*sequences.size()); // everybody is alive
	while (alive > 1) {
		// collapse i and j
#ifdef DEBUG
		cout << "collapse " << max_i << " and " << max_j << endl;
#endif
		profiles[max_i] = profiles[max_i] + profiles[max_j];
		corpse[max_j] = true; // killing max_j
		alive--;
		// updating distances
		// ...in max_i column
		for (unsigned int k = max_i + 1; k < sequences.size(); k++) {
			if (!corpse[k]) {
				if (k > max_j)
					matrix[k*sequences.size()+max_i] += matrix[k*sequences.size()+max_j];
				else 
					matrix[k*sequences.size()+max_i] += matrix[max_j*sequences.size()+k];
				matrix[k*sequences.size()+max_i] /= 2;
			}	
		}
		// ...in max_i row
		for (unsigned int k = 0; k < max_i; k++) {
			if (!corpse[k]) {
				if (k > max_j)
					matrix[max_i*sequences.size()+k] += matrix[k*sequences.size()+max_j];
				else 
					matrix[max_i*sequences.size()+k] += matrix[max_j*sequences.size()+k];
				matrix[max_i*sequences.size()+k] /= 2;
			}
		}
		
#ifdef DEBUG
	for (unsigned int i = 0; i < sequences.size(); i++) {
		if (corpse[i]) continue;
		for (unsigned int j = 0; j < i; j++) {
			if (corpse[j]) continue;
			cout << matrix[i*sequences.size() + j] << '\t';
		}
		cout << '.' << endl;
	}
	cout << endl;
#endif
		
		// searching new max
		max = -MAXINT;
		for (unsigned int i = 0; i < sequences.size(); i++) {
			if (corpse[i]) continue;
			for (unsigned int j = 0; j < i; j++) {
				if (corpse[j]) continue;
				if (matrix[i*sequences.size() + j] > max) {
					max = matrix[i*sequences.size() + j];
					max_i = i;
					max_j = j;
				}
			}
		}
	}
	delete[] matrix;
	delete[] profiles;
	// comment next line if you want to use UPGMA again for other sequences
	// and call UPGMAfree in main program in the end
	UPGMAfree();
}

inline void UPGMAfree() { 
	if (best_mvt) delete[] best_mvt;
	if (scores) delete[] scores;
	best_mvt = NULL;
	scores = NULL;
	dim = 0;
}

Profile& Profile :: operator + (Profile& another) {
	// zero check
	if (another.sequences.size() == 0)
		return *this;
	if (sequences.size() == 0) {
		sequences = another.sequences;
		return *this;
	}
	// reallocate memory (if necessary)
	int new_size = 1;
	int dim1 = sequences[0]->Length()+1, dim2 = another.sequences[0]->Length()+1;
	while (dim1 > new_size || dim2 > new_size) {
		new_size *= 2;
	}
	if (new_size > dim) {
		UPGMAfree();
		dim = new_size;
		best_mvt = new unsigned int[dim*dim];
		scores = new float[dim*dim];
	}
	// initialization
	memset(best_mvt, 0xFF, sizeof(unsigned int)*dim*dim);
	memset(scores, 0, sizeof(float)*dim*dim);
	// calculating new profile
	//  * fill score matrix
#ifdef DEBUG
	cout << "seq1 size: " << dim1 << endl;
	cout << "seq2 size: " << dim2 << endl;
	cout << "allocated matrix: " << dim << 'x' << dim << endl;
#endif
	for (int i = 1; i < dim1; i++) {
		for (int j = 1; j < dim2; j++) {
			/* 25 possible moves
				NT align (x10)
				================
				X | XX | -X | XX
				- | -X | XX | XX
				================
				- | XX | X- | X
				X | X- | XX | X
				================
				  | XX | -- |
				  | -- | XX |
			*/
			float score = ColumnNTscore(another, i-1, j-1) + scores[(i-1)*dim+j-1];
			int way = 15;
			if (j - 1 >= 0) {
				if (best_mvt[i*dim+j-1] < 8 // gap extension
					&& score < scores[i*dim+j-1] + parameters[1]) {
					score = scores[i*dim+j-1] + parameters[1];
					way = 1;
				} else if (best_mvt[i*dim+j-1] > 7 // gap open
					&& score < scores[i*dim+j-1] + parameters[0]) {
					score = scores[i*dim+j-1] + parameters[0];
					way = 1;
				}
			}
			if (i - 1 >= 0) {
				if (best_mvt[(i-1)*dim+j] > 7 && best_mvt[(i-1)*dim+j] < 15 
					&& score < scores[(i-1)*dim+j] + parameters[1]) {
					score = scores[(i-1)*dim+j] + parameters[1];
					way = 8;
				} else if ((best_mvt[(i-1)*dim+j] > 14 || best_mvt[(i-1)*dim+j] < 8)
					&& score < scores[(i-1)*dim+j] + parameters[0]) {
					score = scores[(i-1)*dim+j] + parameters[0];
					way = 8;
				}
			}
			if (i - 2 >= 0 && j - 2 >= 0) {
				if (score < ColumnNTscore(another, i-1, j-1) 
					+ ColumnNTscore(another, i-2, j-2)) {
					score = ColumnNTscore(another, i-1, j-1) 
						+ ColumnNTscore(another, i-2, j-2);
					way = 16;
				}
			}
			if (i - 2 >= 0 && j - 1 >= 0) {
				if (best_mvt[(i-2)*dim+j-1] < 8
					&& score < scores[(i-2)*dim+j-1] + parameters[1] 
					+ ColumnNTscore(another, i-1, j-1)) {
					score = scores[(i-2)*dim+j-1] + parameters[1] 
						+ ColumnNTscore(another, i-1, j-1);
					way = 17;
				} else if (best_mvt[(i-2)*dim+j-1] > 7
					&& score < scores[(i-2)*dim+j-1] + parameters[0] 
					+ ColumnNTscore(another, i-1, j-1)) {
					score = scores[(i-2)*dim+j-1] + parameters[0] 
						+ ColumnNTscore(another, i-1, j-1);
					way = 17;
				} 
				if (score < scores[(i-2)*dim+j-1] + ColumnNTscore(another, i-2, j-1) 
					+ parameters[0]) {
					score = scores[(i-2)*dim+j-1] + ColumnNTscore(another, i-2, j-1) 
						+ parameters[0];
					way = 2;
				}
			}
			if (i - 1 >= 0 && j - 2 >= 0) {
				if (best_mvt[(i-1)*dim+j-2] > 7 && best_mvt[(i-1)*dim+j-2] < 15
					&& score < scores[(i-1)*dim+j-2] + parameters[1] 
					+ ColumnNTscore(another, i-1, j-1)) {
					score = scores[(i-1)*dim+j-2] + parameters[1] 
						+ ColumnNTscore(another, i-1, j-1);
					way = 18;
				} else if ((best_mvt[(i-1)*dim+j-2] > 14 || best_mvt[(i-1)*dim+j-2] < 8)
					&& score < scores[(i-1)*dim+j-2] + parameters[0] 
					+ ColumnNTscore(another, i-1, j-1)) {
					score = scores[(i-1)*dim+j-2] + parameters[0] 
						+ ColumnNTscore(another, i-1, j-1);
					way = 18;
				} 
				if (score < scores[(i-1)*dim+j-2] + ColumnNTscore(another, i-1, j-2) 
					+ parameters[0]) {
					score = scores[(i-1)*dim+j-2] + ColumnNTscore(another, i-1, j-2) 
						+ parameters[0];
					way = 9;
				}
			}
			if (i - 2 >= 0) {
				if (best_mvt[(i-2)*dim+j] < 8 
					&& score < scores[(i-2)*dim+j] + parameters[1] * 2) {
					score = scores[(i-2)*dim+j] + parameters[1] * 2;
					way = 3;
				} else if (best_mvt[(i-2)*dim+j] > 7 
					&& score < scores[(i-2)*dim+j] + parameters[0] + parameters[1]) {
					score = scores[(i-2)*dim+j] + parameters[0] + parameters[1];
					way = 3;
				}
			}
			if (j - 2 >= 0) {
				if (best_mvt[i*dim+j-2] > 7 && best_mvt[i*dim+j-2] < 15
					&& score < scores[i*dim+j-2] + parameters[1] * 2) {
					score = scores[i*dim+j-2] + parameters[1] * 2;
					way = 10;
				} else if ((best_mvt[i*dim+j-2] > 14 || best_mvt[i*dim+j-2] < 8)
					&& score < scores[i*dim+j-2] + parameters[0] + parameters[1]) {
					score = scores[i*dim+j-2] + parameters[0] + parameters[1];
					way = 10;
				} 
			}
			score += parameters[2] * 2;
			// AA align
			if (i - 3 >= 0 && j - 3 >= 0) { 
#ifdef DEBUG
cout << i << ' ' << j << ' ' << score << ' ' << way << ' ' << ColumnAAscore(another, i-3,j-3) << ' ' << ColumnNTscore(another, i-3,j-3) << ' ' << ColumnNTscore(another, i-2,j-2) << ' ' << ColumnNTscore(another, i-1,j-1) << endl;
#endif
				if (score < ColumnAAscore(another, i-3,j-3) 
					+ ColumnNTscore(another, i-3,j-3) + ColumnNTscore(another, i-2,j-2) 
					+ ColumnNTscore(another, i-1, j-1) + scores[(i-3)*dim+j-3]) {
					score = ColumnAAscore(another, i-3,j-3) 
					+ ColumnNTscore(another, i-3,j-3) + ColumnNTscore(another, i-2,j-2) 
					+ ColumnNTscore(another, i-1, j-1)  + scores[(i-3)*dim+j-3];
					way = 19;
				}
			}
			if (i - 3 >= 0 && j - 2 >= 0) {
				if (best_mvt[(i-3)*dim+j-2] < 8) {
					if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(another, i-2, j-2) 
						+ ColumnNTscore(another, i-1, j-1) + parameters[1]) {
						score = scores[(i-3)*dim+j-2] + ColumnNTscore(another, i-2, j-2) 
						+ ColumnNTscore(another, i-1, j-1) + parameters[1];
						way = 20;
					}
				} else {
					if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(another, i-2, j-2) 
						+ ColumnNTscore(another, i-1, j-1) + parameters[0]) {
						score = scores[(i-3)*dim+j-2] + ColumnNTscore(another, i-2, j-2) 
						+ ColumnNTscore(another, i-1, j-1) + parameters[0];
						way = 20;
					}
				}
				if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-1, j-1) + parameters[0]) {
					score = scores[(i-3)*dim+j-2] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-1, j-1) + parameters[0];
					way = 22;
				}
				if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-2, j-1) + parameters[0]) {
					score = scores[(i-3)*dim+j-2] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-2, j-1) + parameters[0];
					way = 4;
				}
			}
			if (i - 3 >= 0 && j - 1 >= 0) {
				int prev_gap = (best_mvt[(i-3)*dim+j-1] < 8) ? 1 : 0;
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-1, j-1)
					+ parameters[prev_gap] + parameters[1]) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-1, j-1)
					+ parameters[prev_gap] + parameters[1];
					way = 23;
				}
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-2, j-1)
					+ parameters[prev_gap] + parameters[0]) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-2, j-1)
					+ parameters[prev_gap] + parameters[0];
					way = 5;
				}
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-3, j-1)
					+ parameters[0] + parameters[1]) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-3, j-1)
					+ parameters[0] + parameters[1];
					way = 6;
				}
			}
			if (i - 3 >= 0) {
				if (best_mvt[(i-3)*dim+j-1] < 8) {
					if (score < scores[(i-3)*dim+j] + parameters[1] * 3) {
						score = scores[(i-3)*dim+j] + parameters[1] * 3;
						way = 7;
					}
				} else {
					if (score < scores[(i-3)*dim+j] + parameters[1] * 2 + parameters[0]) {
						score = scores[(i-3)*dim+j] + parameters[1] * 2 + parameters[0];
						way = 7;
					}
				}
			}
			//
			if (i - 2 >= 0 && j - 3 >= 0) {
				if (best_mvt[(i-2)*dim+j-3] > 7 && best_mvt[(i-2)*dim+j-3] < 15) {
					if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-2) 
						+ ColumnNTscore(another,i-1, j-1) + parameters[1]) {
						score = scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-2) 
						+ ColumnNTscore(another,i-1, j-1) + parameters[1];
						way = 21;
					}
				} else {
					if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-2) 
						+ ColumnNTscore(another,i-1, j-1) + parameters[0]) {
						score = scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-2) 
						+ ColumnNTscore(another,i-1, j-1) + parameters[0];
						way = 21;
					}
				}
				if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-3) 
					+ ColumnNTscore(another,i-1, j-1) + parameters[0]) {
					score = scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-1, j-1) + parameters[0];
					way = 24;
				}
				if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-3) 
					+ ColumnNTscore(another,i-1, j-2) + parameters[0]) {
					score = scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-2, j-1) + parameters[0];
					way = 11;
				}
			}
			if (i - 1 >= 0 && j - 3 >= 0) {
				int prev_gap = (best_mvt[(i-1)*dim+j-3] > 7 
					&& best_mvt[(i-1)*dim+j-3] < 15) ? 1 : 0;
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-1)
					+ parameters[prev_gap] + parameters[1]) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-1)
					+ parameters[prev_gap] + parameters[1];
					way = 25;
				}
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-2)
					+ parameters[prev_gap] + parameters[0]) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-2)
					+ parameters[prev_gap] + parameters[0];
					way = 12;
				}
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-3)
					+ parameters[0] + parameters[1]) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-3)
					+ parameters[0] + parameters[1];
					way = 13;
				}
			}
			if (j - 3 >= 0) {
				if (best_mvt[i*dim+j-3] > 7 && best_mvt[i*dim+j-3] < 15) {
					if (score < scores[i*dim+j-3] + parameters[1] * 3) {
						score = scores[i*dim+j-3] + parameters[1] * 3;
						way = 14;
					}
				} else {
					if (score < scores[i*dim+j-3] + parameters[1] * 2 + parameters[0]) {
						score = scores[i*dim+j-3] + parameters[1] * 2 + parameters[0];
						way = 14;
					}
				}
			}
			// score & way contain best option
			scores[i*dim+j] = score;
			best_mvt[i*dim+j] = way;
		}
	}
	//  * update sequences
	// searching the best score
	int result_score = scores[(dim1-1)*dim+dim2-1];
	int i = dim1-1, j = dim2-1;
	// ...search in last line
	for (int index = 0; index < dim2-1; index++) 
		if (scores[(dim1-1)*dim+index] > result_score) {
			j = index;
			result_score = scores[(dim-1)*dim2+index];
		}
	// ...search in last column
	for (int index = 0; index < dim1-1; index++) 
		if (scores[index*dim + dim2-1] > result_score) {
			i = index; 
			j = dim2-1;
			result_score = scores[index*dim + dim2-1];
		}
	// i, j - point of best score
	// updating sequences
	this->InsertGap(dim1-1, dim2-1-j); // insert gaps in Profile 1
	another.InsertGap(dim2-1, dim1-1-i); // insert gaps in Profile 2
	
#ifdef DEBUG
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			cout << best_mvt[i*dim+j] << '\t';
		}
		cout << endl;
	}
	cout << endl << endl;
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			cout << scores[i*dim+j] << '\t';
		}
		cout << endl;
	}
	cout << endl;
#endif

	
	
	while (i && j) {
		switch (best_mvt[i*dim + j]) {
			case 1:
				another.InsertGap(j);
				i--;
				break;
			case 2:
				another.InsertGap(j);
				i -= 2;
				j--;
				break;
			case 3:
				another.InsertGap(j, 2);
				i -= 2;
				break;
			case 4:
				another.InsertGap(j);
				i -= 3;
				j -= 2;
				break;
			case 5:
				another.InsertGap(j);
				j--;
				another.InsertGap(j);
				i -= 3;
				break;
			case 6:
				another.InsertGap(j, 2);
				i -= 3;
				j--;
				break;
			case 7:
				another.InsertGap(j, 3);
				i -= 3;
				break;
			case 8:
				InsertGap(i);
				j--;
				break;
			case 9:
				InsertGap(i);
				i--;
				j -= 2;
				break;
			case 10:
				InsertGap(i, 2);
				j -= 2;
				break;
			case 11:
				InsertGap(i);
				i -= 2;
				j -= 3;
				break;
			case 12:
				InsertGap(i);
				i--;
				InsertGap(i);
				j -= 3;
				break;
			case 13:
				InsertGap(i, 2);
				i--;
				j -= 3;
				break;
			case 14:
				InsertGap(i, 3);
				j -= 3;
				break;
			case 15:
				i--;
				j--;
				break;
			case 16:
				i -= 2;
				j -= 2;
				break;
			case 17:
				i -= 2;
				j--;
				another.InsertGap(j);
				break;
			case 18:
				i--;
				InsertGap(i);
				j -= 2;
				break;
			case 19:
				i -= 3;
				j -= 3;
				break;
			case 20:
				i -= 3;
				j -= 2;
				another.InsertGap(j);
				break;
			case 21:
				i -= 2;
				InsertGap(i);
				j -= 3;
				break;
			case 22:
				i -= 3;
				j--;
				another.InsertGap(j);
				j--;
				break;
			case 23:
				i -= 3;
				j--;
				another.InsertGap(j, 2);
				break;
			case 24:
				i--;
				InsertGap(i);
				i--;
				j -= 3;
				break;
			case 25:
				i--;
				InsertGap(i, 2);
				j -= 3;
				break;
			default:
				// some problems
				// return NULL; // DEBUG
				return *this;
		}
	}
	InsertGap(0, j);
	another.InsertGap(0, i);
	return *this;
}

float Profile :: ColumnNTscore(Profile& another, int index1, int index2) {
	float score = 0;
	float denominator = sequences.size() + another.sequences.size();
	CalcFrequenciesNT(index1); 
	another.CalcFrequenciesNT(index2);
	for (unsigned int i = 0; i < sequences.size(); i++) {
		unsigned char char1 = sequences[i]->nt_seq[index1];
		if (char1 == '-') continue;
		for (unsigned int j = 0; j < another.sequences.size(); j++) {
			unsigned char char2 = another.sequences[j]->nt_seq[index2];
			if (char2 == '-') continue;
			float numerator = frequency[char1]+another.frequency[char2];
			score += nt_score_matrix[char1*128+char2]*(numerator/denominator);
		}
	}
	return score;
}

float Profile :: ColumnAAscore(Profile& another, int index1, int index2) {
	float score = 0;
	float denominator = sequences.size() + another.sequences.size();
	CalcFrequenciesAA(index1); 
	another.CalcFrequenciesAA(index2);
	for (unsigned int i = 0; i < sequences.size(); i++) {
		unsigned char char1 = sequences[i]->TranslateNTtoAA(index1);
		if (char1 == '-' || char1 == '!' || char1 == '*') continue;
		for (unsigned int j = 0; j < another.sequences.size(); j++) {
			unsigned char char2 = another.sequences[j]->TranslateNTtoAA(index2);
			if (char2 == '-' || char2 == '!' || char2 == '*') continue;
			float numerator = frequency[char1]+another.frequency[char2];
			score += aa_score_matrix[char1*128+char2]*(numerator/denominator);
		}
	}
	float addition1 = 2 * (frequency['!'] + another.frequency['!']) / denominator;
	addition1 *= parameters[2]; // gap frame
	float addition2 = 2 * (frequency['*'] + another.frequency['*']) / denominator;
	addition2 *= parameters[3]; // stop cost
	return score + addition1 + addition2;
}

void Profile :: CalcFrequenciesAA(int position) {
	memset(frequency, 0, sizeof(int)*128);
	for_each(sequences.begin(), sequences.end(), 
		[&](BioSeq* s) { 
			frequency[(unsigned char)s->TranslateNTtoAA(position)]++; 
		} );
}

void Profile :: CalcFrequenciesNT(int position) {
	memset(frequency, 0, sizeof(int)*128);
	for_each(sequences.begin(), sequences.end(), 
		[&](BioSeq* s) { 
			frequency[(unsigned char)s->nt_seq[position]]++; 
		} );
}
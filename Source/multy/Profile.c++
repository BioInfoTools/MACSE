#include <string.h>

#include "Profile.h"
#include "PairwiseAlign.h"

#define MAXINT 0x7FFFFFFF

// global tables for merge profiles
// memory allocation in Profile::operator +
// free in UPGMAfree
int* best_mvt = NULL;
int* scores = NULL;
int dim = 0;

void UPGMA(std::vector<BioSeq*>& sequences, PairwiseAlign& aligner) {
	int* matrix = new int[sequences.size() * sequences.size()];
	Profile* profiles = new Profile[sequences.size()];
	int max = -MAXINT;
	unsigned int max_i, max_j;
	for (unsigned int i = 0; i < sequences.size(); i++) {
		// at first each sequence is a single cluster
		profiles[i].sequences.push_back(sequences[i]);
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
	int alive = sequences.size();
	bool* corpse = new bool[sequences.size()];
	memset(corpse, false, sizeof(bool)*sequences.size()); // everybody is alive
	while (alive > 1) {
		// collapse i and j
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
}

inline void UPGMAfree() { 
	if (best_mvt) delete[] best_mvt;
	if (scores) delete[] scores;
	dim = 0;
}

Profile& Profile :: operator + (const Profile& another) {
	// zero check
	if (another.sequences.size() == 0)
		return *this;
	if (sequences.size() == 0) {
		sequences = another.sequences;
		return *this;
	}
	// reallocate memory (if necessary)
	int new_size = 1;
	int dim1 = sequences[0]->Length(), dim2 = another.sequences[0]->Length();
	while (dim1 + 1 > new_size || dim2 + 1 > new_size) {
		new_size *= 2;
	}
	if (new_size > dim) {
		UPGMAfree();
		best_mvt = new int[dim*dim];
		scores = new int[dim*dim];
		dim = new_size;
	}
	// initialization
	memset(best_mvt, 0, sizeof(int)*dim*dim);
	memset(scores, 0, sizeof(int)*dim*dim);
	// calculating new profile
	//  * fill score matrix
	for (int i = 1; i < dim1; i++) {
		for (int j = 1; j < dim2; j++) {
			// 27 possible moves
			
		}
	}
	//  * update sequences
	
	return *this;
}

int Profile :: ColumnNTscore(const Profile& another, int index) {
	int score = 0;
	for (unsigned int i = 0; i < sequences.size(); i++) {
		char char1 = sequences[i]->nt_seq[index];
		if (char1 == '-') continue;
		for (unsigned int j = 0; j < another.sequences.size(); j++) {
			char char2 = another.sequences[j]->nt_seq[index];
			if (char2 == '-') continue;
			score += nt_score_matrix[char1*128+char2];
		}
	}
	return score;
}

int Profile :: ColumnAAscore(const Profile& another, int index) {
	int score = 0;
	for (unsigned int i = 0; i < sequences.size(); i++) {
		char char1 = sequences[i]->TranslateNTtoAA(index);
		if (char1 == '-' || char1 == '!') continue;
		for (unsigned int j = 0; j < another.sequences.size(); j++) {
			char char2 = another.sequences[j]->TranslateNTtoAA(index);
			if (char2 == '-' || char2 == '!') continue;
			score += aa_score_matrix[char1*128+char2];
		}
	}
	return score;
}
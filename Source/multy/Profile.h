#ifndef _Profile_H_
#define _Profile_H_

#include <vector>
#include <algorithm>

#include "BioSeq.h"
#include "PairwiseAlign.h"

struct Profile {
	// data
	std::vector<BioSeq*> sequences;
	const int* nt_score_matrix = NULL; // int array [128*128]
	const int* aa_score_matrix = NULL; // int array [128*128]
	// constructor
	Profile() { };
	// destructor
	~Profile() {
		// everything will be deleted in main program
		/*
		for_each(sequences.begin(), sequences.end(), 
			[] (BioSeq* seq) { delete seq; } );
		delete[] nt_score_matrix;
		delete[] aa_score_matrix;
		*/
	};
	// alignment of two profiles
	Profile& operator + (const Profile& another);
	// get score in column
	int ColumnNTscore(const Profile& another, int index);
	int ColumnAAscore(const Profile& another, int index);
};

// UPGMA function change input sequences!
void UPGMA(std::vector<BioSeq*>& sequences, PairwiseAlign& aligner);
void UPGMAfree();

#endif
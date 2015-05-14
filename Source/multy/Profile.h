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
	int frequency[128];
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
	Profile& operator + (Profile& another);
	// get score in column
	float ColumnNTscore(Profile& another, int index1, int index2);
	float ColumnAAscore(Profile& another, int index1, int index2);
	// fill an array of frequency
	void CalcFrequenciesNT(int position);
	void CalcFrequenciesAA(int position);
	// insert gaps in sequences
	void InsertGap(int pos, int count) {
		for_each(sequences.begin(), sequences.end(), [pos, count] (BioSeq* s) {
			if (count) s->nt_seq.insert(pos, count, '-');
		});
	}
	void InsertGap(int pos) {
		for_each(sequences.begin(), sequences.end(), [pos] (BioSeq* s) {
			s->nt_seq.insert(pos, 1, '-');
		});
	}
};

// UPGMA function change input sequences!
void UPGMA(std::vector<BioSeq*>& sequences, PairwiseAlign& aligner);
void UPGMAfree();

#endif
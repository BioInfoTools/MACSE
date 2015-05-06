#ifndef _Profile_H_
#define _Profile_H_

#include <vector>
#include <algorithm>

#include "BioSeq.h"
#include "PairwiseAlign.h"

struct Profile {
	// data
	std::vector<BioSeq*> sequences;
	// constructor
	Profile() { };
	// destructor
	~Profile() {
		/*
		for_each(sequences.begin(), sequences.end(), 
			[] (BioSeq* seq) { delete seq; } );
		*/
		// each sequence will be deleted in main programm!
	};
	// alignment of two profiles
	Profile& operator + (const Profile& another);
};

// UPGMA function change input sequences!
void UPGMA(std::vector<BioSeq*>& sequences, PairwiseAlign& aligner);

#endif
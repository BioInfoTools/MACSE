#include <algorithm> 
#include <stdlib.h>

#include "BioSeq.h"
#include "Reader.h"
#include "Profile.h"
#include "PairwiseAlign.h"

#define HELP "./multy <Seq.fasta> [<NTsubs> <AAsubs> <GapOpen> <GapExt> <GapFrame> <StopCost>]" 
#define WARN "WARNING: wrong arguments count; using the default values"

using std::cout;
using std::endl;

int main(int argc, char** argv) {
	if (argc > 1) {
		// input fasta file to align
		char* input_fasta = argv[1];
		// default values
		std::string NTsubs = "NUC-45", AAsubs = "BLOSUM62";
		int gap_open = -2, gap_extension = -1, gap_frame = -2, stop_cost = -3;
		// extra parameters
		if (argc == 8) {
			NTsubs = argv[2];
			AAsubs = argv[3];
			gap_open = atoi(argv[4]);
			gap_extension = atoi(argv[5]);
			gap_frame = atoi(argv[6]);
			stop_cost = atoi(argv[7]);
		} else if (argc != 2) cout << WARN << endl;
		//===============
		cout << "Reading sequences..." << endl;
		std::vector<BioSeq*> data;
		ReadFastaFile(input_fasta, data);
		cout << data.size() << " sequences were obtained" << endl; 
		cout << "Input parameters:" << endl;
		cout << "NT substitution matrix\t" << NTsubs << endl;
		cout << "AA substitution matrix\t" << AAsubs << endl;
		cout << "Gap open cost\t\t" << gap_open << endl;
		cout << "Gap extension cost\t" << gap_extension << endl;
		cout << "Gap frame cost\t\t" << gap_frame << endl;
		cout << "Stop codon cost\t\t" << stop_cost << endl;
		// creating aligner tool
		PairwiseAlign aligner(NTsubs.c_str(), AAsubs.c_str(), stop_cost, gap_open, 
													gap_extension, gap_frame);
		// calculating MSA
		UPGMA(data, aligner); // this call will change data!
		// print results
		cout << endl;
		for_each(data.begin(), data.end(), [] (BioSeq* seq) {seq->PrintNT(cout);} );
		cout << endl;
		for_each(data.begin(), data.end(), [] (BioSeq* seq) {seq->PrintAA(cout);} );
		// free memory
		for_each(data.begin(), data.end(), [] (BioSeq* seq) { delete seq; } );
	} else {
		std::cout << "ERROR: you must type an input file name" << endl;
		std::cout << HELP << std::endl;
	}
	return 0;
}
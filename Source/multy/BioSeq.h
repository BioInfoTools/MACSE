#ifndef _BioSeq_H_
#define _BioSeq_H_

#include <string>
#include <iostream>

struct BioSeq {
	// data
	std::string name;
	std::string nt_seq;
	int first_nongap = 0;
	int last_nongap  = 0;
	// constructor
	BioSeq(std::string n, std::string s): name(n), nt_seq(s) { }
	// destructor
	~BioSeq() { }
	// accessors
	int Length() { return nt_seq.length(); }
	char operator [] (int index) { return nt_seq[index]; }
	// translater
	std::string GetAAseq();
	char TranslateNTtoAA(int index) const;
	// printer
	void PrintAA(std::ostream& out);
	void PrintNT(std::ostream& out) { 
		out << name   << std::endl;
		out << nt_seq << std::endl; 
	}
};

#endif
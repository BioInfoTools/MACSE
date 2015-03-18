#include <string>

void BLOSUM62(int* matrix);
void BLOSUM50(int* matrix);
void NTABLE(int* matrix);

typedef std::pair<std::string, std::string> string_tuple;

class AlignMethod {
	protected:
		std::string seq1, seq2;
		std::string align1, align2;
		int score_matrix[128*128];
		int result_score;
	public:
		virtual ~AlignMethod() { }
		virtual string_tuple Align()=0;
		virtual void NewStrings(std::string s1, std::string s2)=0;
		string_tuple GetAlign() { return std::make_pair(align1, align2); }
		int GetScore() { return result_score; }
};

class NeedlmanWunsch: public AlignMethod {
	private:
        int gap_open, gap_extension;
		int* F = NULL;
		int* W = NULL;
		int n, m;
	public:
        NeedlmanWunsch(std::string s1, std::string s2,
                       std::string score_matrix,
                       int gap_open, int gap_extension);
        virtual ~NeedlmanWunsch() { if (F) delete [] F; if (W) delete [] W; }
        //====================
        string_tuple Align();
        void NewStrings(std::string s1, std::string s2);
        //====================
		void NewScoreMatrix(std::string score_matrix);
		void ChangeGapOpen(int new_value) { gap_open = new_value; }
		void ChangeGapExtension(int new_value) { gap_extension = new_value; }
};

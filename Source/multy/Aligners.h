#ifndef ALIGNERS_H
#define ALIGNERS_H

#include <string>

typedef std::pair<std::string, std::string> string_tuple;

class AlignMethod {
	protected:
        //общие элементы методов выравнивания
		std::string seq1, seq2;
		std::string align1, align2;
		int gap_open, gap_extension;
		int result_score;
		//замена матрицы совпадений
		void NewScoreMatrix(std::string file_name, int* score_matrix);
	public:
		virtual ~AlignMethod() { }
		//метод выравнивания
		virtual string_tuple Align()=0;
		//аксессоры
		int GetScore() { return result_score; }
		string_tuple GetAlign() { return std::make_pair(align1, align2); }
		//модификаторы
		void ChangeGapOpen(int new_value) { gap_open = new_value; }
		void ChangeGapExtension(int new_value) { gap_extension = new_value; }
};

class NeedlmanWunsch: public AlignMethod {
	private:
        //матрица замен
        int score_matrix[128*128];
        //матрицы для ДП
		int* F = NULL;  //матрица очков
		int* W = NULL;  //матрица обратного хода
		int n, m; // размеры строк +1
	public:
        NeedlmanWunsch(std::string s1, std::string s2,
                       std::string score_matrix,
                       int gap_open, int gap_extension);
        virtual ~NeedlmanWunsch() { if (F) delete [] F; if (W) delete [] W; }
        //====================
        string_tuple Align();
        //====================
		void ChangeScoreMatrix(std::string file_name) {
            NewScoreMatrix(file_name, score_matrix);
        }
        void ChangeStrings(std::string s1, std::string s2);
};

class Hein: public AlignMethod {
	private:
        int n, m, gap_frame;
		int* F = NULL;
		int* W = NULL;
        std::string AAseq1, AAseq2;
        int nt_score_matrix[128*128];
        int aa_score_matrix[128*128];
        std::string TranslateNTtoAA(std::string& s);
	public:
        Hein(std::string s1, std::string s2,
              std::string nt_score_matrix,
              std::string aa_score_matrix,
              int gap_open, int gap_extension, int frame_gap);
        virtual ~Hein() { }
        //====================
        string_tuple Align();
        //====================
		void ChangeNTscoreMatrix(std::string file_name) {
            NewScoreMatrix(file_name, nt_score_matrix);
		}
		void ChangeAAscoreMatrix(std::string file_name) {
            NewScoreMatrix(file_name, aa_score_matrix);
		}
		void ChangeStrings(std::string s1, std::string s2);
};

#endif

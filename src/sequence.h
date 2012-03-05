#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>

using namespace std;

class Sequence{
private:
	string id;
	string comment;
	string name;
	string seq;
	string al_seq;
	int sqlite_id;
	bool aligned;
	string reverse(string );

public:
	Sequence();
	Sequence(string,string,bool);
	Sequence(string,string);
	bool is_aligned();
	string get_sequence();
	string get_id();//this will be printed for the genbank
	string get_name();
	string get_comment();
	void set_sequence(string seq);
	void set_id(string id);
	void set_name(string name);
	void set_comment(string comment);
	void set_aligned(bool al);
	string reverse_complement();
	void perm_reverse_complement();
	void set_sqlite_id(int iid);
	int get_sqlite_id();
	string get_aligned_seq();
	void set_aligned_seq(string inseq);
};
#endif

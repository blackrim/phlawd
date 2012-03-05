#ifndef _GENEDB_H
#define _GENEDB_H

#include <string>
#include <vector>

#include "utils.h"
#include "DBSeq.h"
#include "sequence.h"

using namespace std;

class GeneDB{
private:
    string name;//database filename
    
public:
    GeneDB();
    GeneDB(string name);
    void add_user_seqs_to_db(vector<Sequence> * user_seqs);
    void add_seqs_to_db(vector<Sequence> * keep_seqs);
    int add_alignment(string filen, vector<Sequence> * dbseqs, vector<Sequence> * userseqs);
    void remove_alignment(int alignid);
    void remove_alignment_by_name(string alignname);
    void add_seq_to_alignment(int alignid,Sequence inseq);
    void get_align_seqs(int alignid, vector<Sequence> & seqs);
    void get_align_seqs_unaligned(int alignid, vector<Sequence> & seqs);
    void get_align_seq_unaligned_fully_initialized(string alignname,vector<Sequence> & seqs);
    void initialize(bool overwrite);
    void update_align_seqs(int alignid,vector<Sequence> & seqs);
    void get_all_sequences(vector<Sequence> & seqs);
    void get_alignment_names(vector<string> & names);
};

#endif

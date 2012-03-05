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
    void add_seqs_to_db(vector<DBSeq> * keep_seqs);
    int add_alignment(string filen, vector<DBSeq> * dbseqs, vector<Sequence> * userseqs);
    void remove_alignment(int alignid);
    void add_seq_to_alignment(int alignid,Sequence inseq);
    void get_align_seqs(int alignid, vector<Sequence> & seqs);
    void initialize(bool overwrite);
    void update_align_seqs(int alignid,vector<Sequence> & seqs);
};

#endif

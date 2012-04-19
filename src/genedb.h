#ifndef _GENEDB_H
#define _GENEDB_H

#include <string>
#include <vector>
#include <map>

#include "utils.h"
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
    int get_alignment_id_by_name(string alignname);
    void toggle_alignment_update(int alignid);
    void remove_alignment_by_name(string alignname);
    void add_seq_to_alignment(int alignid,Sequence inseq);
    void get_align_seqs(int alignid, vector<Sequence> & seqs);
    void get_profile_align_seqs(int alignid, vector<Sequence> & seqs);
    void get_align_seqs_unaligned(int alignid, vector<Sequence> & seqs);
    void get_align_seq_unaligned_fully_initialized(string alignname,vector<Sequence> & seqs);
    void initialize(bool overwrite);
    void update_align_seqs(int alignid,vector<Sequence> & seqs);
    void update_profile_align_seqs(int alignid, vector<Sequence> & seqs);
    void get_all_sequences(vector<Sequence> & seqs);
    void get_first_profile_alignments(vector<string> & names);
    void get_alignment_names(vector<string> & names);
    void get_alignment_nums(vector<int> & nums);
    void get_profile_alignment_nums(vector<int>&);
    void remove_profile_alignments();
    void copy_alignments_to_first_profiles(map<int,string> & profile_id_name_map);
    void copy_alignments_to_first_profiles_updated(map<int, string> & profile_id_name_map,vector<int>& updatedprofsnums);
    int get_deepest_profile_for_alignment(int alignid);
    int add_profile_alignment(int child1, int child2);
    void add_sequences_for_profile_alignment(int profilealignid,vector<Sequence> & seqs);
    void write_profile_alignment_to_file(int alignid,string filename);
    void write_profile_alignment_with_names_to_file(int alignid, string filename,bool ncbi);
    void get_updated_profs_names_delete_old(vector<string> & updatedprofs,vector<int> & updatedprofsnums,vector<int> & notupdatedprofnums);
    void toggle_updated_all_off();
};

#endif

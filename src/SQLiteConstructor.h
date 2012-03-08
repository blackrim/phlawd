/*
 * PHLAWD: Phylogeny assembly with databases
 * Copyright (C) 2010  Stephen A. Smith
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

/*
 *  SQLiteConstructor.h
 */

#ifndef _SQLITECONSTRUCTOR_H_
#define _SQLITECONSTRUCTOR_H_

#include <string>
#include <vector>
#include <map>

using namespace std;

#include "libsqlitewrapped.h"

#include "sequence.h"

#include "tree.h"
#include "node.h"
#include "genedb.h"

class SQLiteConstructor {
private:
    ofstream logfile;
    fstream gifile;
    fstream ufafile;
    string clade_name;
    vector <string> search;
    string gene_name;
    string gene_db_name;
	string genefoldername;
    GeneDB gene_db;
    double mad_cutoff;
    double coverage;
    double identity;
    string db;
    int port;
    string listfilename;
    string excludefilename;
    string exclude_gi_filename;
    string include_gi_filename;
    bool onlynamesfromfile;
    bool excludenamesfromfile;
    bool includegifromfile;
    bool excludegifromfile;
    bool containshigher;
    bool containswild;
    bool useITS;
    int numthreads;
    bool automated;
    bool updateDB;
    bool updateFILE;
    string updatef;
    bool ncbi_saturation;
    bool usertree;
    string usertreefile;
    Tree * userguidetree;
    bool userfasta;
    string userfastafile;
    bool userskipsearchdb;
    bool skipdbcheck;
    int main_left; // left value for clade id
    int main_right; // right value for clade id
    map<Sequence*,Node*> user_fasta_node_map; 
    vector<Sequence> * user_seqs;
    vector<Sequence> * known_seqs;
    vector<Sequence> use_only_names_from_file(vector<Sequence>& seqs);
    Sequence add_higher_taxa(string taxon_id,vector<Sequence>& seqs);
    vector<Sequence> exclude_names_from_file(vector<Sequence>& seqs);
    vector<Sequence> exclude_gis_from_file(vector<Sequence>& seqs);
    vector<Sequence> include_gis_from_file(vector<Sequence>& seqs);
    void first_seq_search_for_gene_left_right(vector<vector<string> >  &);
    vector<Sequence> first_get_seqs_for_name_use_left_right(int name_id, vector<vector<string> > & results);
    void get_same_seqs_openmp_SWPS3(vector<Sequence> & seqs,  vector<Sequence> * keep_seqs);
    double get_usertree_keepseq_overlap(vector<Sequence> * keep_seqs);
    void remove_duplicates_SWPS3(vector<Sequence> * keep_seqs);
    void reduce_genomes(vector<Sequence> * keep_seqs);
    void get_seqs_for_names(string name_id, vector<Sequence> * seqs, vector<Sequence> * temp_seqs);
    void get_seqs_for_names_user(string name_id, vector<Sequence> * seqs);
    void get_seqs_for_nodes(Node * node, vector<Sequence> * seqs, vector<Sequence> * temp_seqs);
    void get_seqs_for_user_nodes(Node * node, vector<Sequence> * seqs);
    vector<string> get_final_children(string name_id);
    vector<string> get_final_children_node(Node * node);
    vector<string> get_final_children_node_hier(Node * node);
    void make_mafft_multiple_alignment(vector<Sequence> * inseqs,vector<Sequence> * inseqs2);
    double calculate_MAD_quicktree();
    double calculate_MAD_quicktree_sample(vector<Sequence> * inseqs, vector<Sequence> * inuserseqs);
    void saturation_tests(vector<string> name_ids, vector<string> names, vector<Sequence> * keep_seqs);
    int get_single_to_group_seq_score(Sequence & inseq,vector<Sequence> & ginseqs);
    void write_gi_numbers(vector<Sequence> *);
    void write_user_numbers();
    void add_seqs_from_db_to_seqs_vector(string alignname,vector<Sequence> * keep_seqs, vector<Sequence> & storedseqs);
    void remove_seq_from_seq_vector(vector<Sequence> * inseqs,string sid);
    void match_aligned_file(vector<Sequence> * temp_seqs, vector<Sequence> * temp_user_seqs);
    void match_an_aligned_seq(Sequence * temp_seq);
    void get_aligned_file(vector<Sequence> * temp_seqs);
public:
    SQLiteConstructor(string cn, vector <string> searchstr, string genen, string genedb,
		      double mad_cut,double cover, double ident, string dbs,
		      string known_seq_filen, bool its, int numt,bool autom,
		      bool updb, string updf);
    ~SQLiteConstructor(){}
    void set_only_names_from_file(string filename, bool containshi, bool containswild);
    void set_exclude_names_from_file(string filename);
    void set_exclude_gi_from_file(string filename);
    void set_include_gi_from_file(string filename);
    void set_user_guide_tree(string filename, bool skipcheckdb);
    void set_user_fasta_file(string filename, bool skipcheckdb);
    void set_user_skip_search();
    int run(); //0 = success, 2 = update>20%
    string get_cladename();
    vector <string> get_search();
    string get_genename();
    string get_genedb();
    double get_madcutoff();
    double get_coverage();
    double get_identity();
    int get_numthreads();
    bool get_updatestatus();
    Tree * get_user_guide_tree_obj();
};


#endif	//_CONSTRUCTOR_H_

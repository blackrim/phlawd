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
 * SQLiteProfiler.h
 */

#ifndef SQLITEPROFILER_H_
#define SQLITEPROFILER_H_

#include "tree.h"
#include "genedb.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

#include "libsqlitewrapped.h"

class SQLiteProfiler{
private:
    string gene_name;
    string gene_db_name;
    GeneDB gene_db;
    vector<string> align_names;
    vector<int> align_nums;
    string cladename;
    string db;
    string profilefoldername;
    bool use_orphan;
    bool automated;
    bool updatedb;
    Tree * userguidetree;
    bool usertree;
    map<int, string> profile_id_name_map;
    vector<int> updatednums;
    vector<int> updatedprofiles;

//functions
    void get_children(string in_id, vector<string> * in_ids, vector<string> * in_keepids);
    vector<string> get_final_children(string id);
    string get_right_one(vector<string> allids,Query & res);
    vector<string> get_left_right_children(string id);
    void create_distances(string clade_name,map<int, map<int,double> > * numlist);
    void create_distances_user_tree(vector<string> names,map<string,string> * numnames
				    ,map<string,string> * namesnum, map<int, map<int,double> > * numlist);
    void get_shortest_distance_with_dicts(vector<int>& nums,map<int, map<int,double> > & numlist, int * shortestnameone, 
					  vector<int> * shortestnametwo);
    void clean_before_profile(int infile);
    int profile(map<int, map<int,double> > numlist);
    void copy_final_file(string filename);
    string get_name_from_tax_id(string taxid);
    void calculate_for_removal(vector<Sequence> * seqs,
			       map<string,double> & allmeans);
    void calculate_for_removal_quicktree(vector<Sequence> * seqs, map<string,double> & allmeans);
    void remove_outliers();
    void rename_final_alignment(int alignid);
    int make_muscle_profile(int profile1,int profile2,int outprofile);
    double get_muscle_spscore(string filename);
    void test_outfile_exists(string filename);
    void match_and_add_profile_alignment_to_db(int profileid);
public:
    SQLiteProfiler(string gn, string gene_dbn,string cn, string dbs, bool autom,bool updb);
    void prelimalign();
    void run();
    void set_user_guide_tree(Tree * tree);
};

#endif /* MQPROFILER_H_ */

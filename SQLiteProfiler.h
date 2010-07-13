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

#include <Phyl/trees>
#include <Phyl/Newick.h>

using namespace bpp;

#include <string>
#include <map>
#include <vector>

using namespace std;

#include "libsqlitewrapped.h"

class SQLiteProfiler{
	private:
		string gene_name;
		vector<string> file_names;
		string cladename;
		string db;
		string profilefoldername;
		bool use_orphan;
		bool automated;
		Tree * read_user_guide_tree(string filen);
		void get_children(string in_id, vector<string> * in_ids, vector<string> * in_keepids);
		vector<string> get_final_children(string id);
		int count_seqs(string file_name);
		string get_right_one(vector<string> allids,Query & res);
		vector<string> get_left_right_children(string id);
		void create_distances(string clade_name, vector<string> names,map<string,string> * numnames,
				map<string,string>* namesnum, vector< vector<double> > * numlist);
		void create_distances(vector<string> names,Tree * tree,map<string,string> * numnames
				,map<string,string> * namesnum, vector< vector<double> > * numlist);
		void get_shortest_distance_with_dicts(vector<string> names,map<string,string> numnames,map<string,string> namesnum,
				vector< vector<double> > numlist, string * shortestnameone, vector<string> * shortestnametwo);
		void clean_before_profile(string infile);
		void profile(map<string,string> numnames,map<string,string> namesnum,
						vector< vector<double> > numlist);
		void copy_final_file(string filename);
		string get_name_from_tax_id(string taxid);
		void calculate_for_removal(SequenceContainer * seqs,
				map<string,double> & allmeans);
		void calculate_for_removal_quicktree(SequenceContainer * seqs,
						map<string,double> & allmeans);
		void remove_outliers();
		void rename_final_alignment(string which);
	public:
		SQLiteProfiler(string gn, string cn, string dbs, bool autom);
		void prelimalign();
		void run();
};

#endif /* MQPROFILER_H_ */

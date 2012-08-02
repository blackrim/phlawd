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
 * SQLiteProfiler.cpp
 */

#include <string.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <unistd.h>

using namespace std;

#include "tree.h"
#include "node.h"
#include "tree_utils.h"
#include "tree_reader.h"
#include "sequence.h"
#include "fasta_util.h"


#include "utils.h"

#include "omp.h" 
#include "libsqlitewrapped.h"

#include "SQLiteProfiler.h"

#define SIXES 666666

template <class T>
inline std::string to_string (const T& t){
    std::stringstream ss;
    ss << t;
    return ss.str();
}

SQLiteProfiler::SQLiteProfiler(string gn, string gene_dbn,string cn, string dbs,bool autom,bool updb): 
                              gene_name(gn),gene_db_name(gene_dbn),
			      use_orphan(false), cladename(cn),db(dbs),automated(autom),updatedb(updb),
			      usertree(false){
    profilefoldername = gene_name+"_TEMPFILES/";
    gene_db= GeneDB(gene_db_name);
}

void SQLiteProfiler::prelimalign() {
    mkdir(profilefoldername.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH);
    bool standard = true;
    updatednums = vector<int> ();
    if(updatedb == true){
	vector<string> samefiles; //files were not updated at all
	bool newfiles=false; // if there are new files from MAD in align stage, kick out and tree as new for now 
	vector<string> originalalnfiles;
	vector<string> curproffiles;
	gene_db.load_orig_alignment_names_into(originalalnfiles);
	gene_db.get_first_profile_alignments(curproffiles);
	if (originalalnfiles.size() != curproffiles.size() || originalalnfiles.size() == 1)
	    newfiles = true;
	else{//double check that there aren't new ones
	    std::sort(originalalnfiles.begin(), originalalnfiles.end());
	    std::sort(curproffiles.begin(), curproffiles.end());
	    std::vector<string> v3;
	    std::set_intersection(originalalnfiles.begin(), originalalnfiles.end(), curproffiles.begin(), curproffiles.end(), std::back_inserter(v3));
	    if (v3.size() != originalalnfiles.size())
		newfiles = true;
	}
	if (newfiles == true){
	    gene_db.remove_profile_alignments();
	    standard = true;
	    updatedb = false;
	}else{//update the profiles
	    standard = false;
	    vector<string> updatedprofs;
	    vector<int> updatedprofsnums;//these will be the profile ids that are updated
	    vector<int> notupdatedprofsnums;//these will not be updated
	    //TODO: start editing here, need the 
	    gene_db.get_updated_profs_names_delete_old(updatedprofs,updatedprofsnums,notupdatedprofsnums);
	    updatedprofiles.clear();
	    for (int i=0;i<notupdatedprofsnums.size();i++){
		int tid = gene_db.get_deepest_profile_for_alignment(notupdatedprofsnums[i]);
		if (tid == -1){
		    updatednums.push_back(notupdatedprofsnums[i]);
		    cout << "id: " << notupdatedprofsnums[i] << endl;
		}
		else if(count(updatedprofiles.begin(),updatedprofiles.end(),tid)==0){
		    updatedprofiles.push_back(tid);
		    cout << "tid: " << tid << endl;
		}
	    }
	}
    }
    if(standard == true){ //not an update run

        // these are properties of the SQLiteProfiler class, accessed by various methods throughout 
        orig_align_names = vector<string>();
        first_profiles_dbids = vector<int>();

        cout << "getting alignments" << endl;
        gene_db.load_orig_alignment_names_into(orig_align_names);
        profile_id_name_map = map<int,string>();
        gene_db.copy_alignments_to_first_profiles(profile_id_name_map);
        gene_db.load_first_profile_ids_into(first_profiles_dbids);

    } else { // updaterun

        // only copy over the ones that are updated
        // only align those files that are updated
        orig_align_names = vector<string>();
        first_profiles_dbids = vector<int>();
        cout << "getting updated alignment files" << endl;
        profile_id_name_map = map<int,string>();

        gene_db.copy_alignments_to_first_profiles_updated(profile_id_name_map,updatednums);
        gene_db.load_orig_alignment_names_into(orig_align_names);
        gene_db.load_first_profile_ids_into(first_profiles_dbids);
    }
}

void SQLiteProfiler::set_user_guide_tree(Tree * tree){
    usertree = true;
    userguidetree = tree;
}

void SQLiteProfiler::run(){
    cout << "starting run" << endl;
    int finalaln;
    if(first_profiles_dbids.size() > 1){
	map<string,string> numnames;
	map<string,string> namesnum;
	map<int, map<int,double> > numlist;
	//TODO: make sure that this works with update , incomplete user guide tree
	if(usertree == true){
	    cout << "user guide tree" << endl;
	    create_distances_user_tree(orig_align_names,&numnames,&namesnum,&numlist);
	}else{//use ncbi tree
	    cout << "ncbi guide tree" << endl;
	    create_distances(cladename,&numlist);
	}
	//start profiling
	cout<<"profiling"<<endl;
	finalaln = profile(numlist);
    }else{
	finalaln = 1;
	clean_dbseqs(1);
    }
    if (updatedb == true)
	gene_db.toggle_updated_all_off();
    cout << "writing final alignment" << endl;
    rename_final_alignment(finalaln);//requires FINAL.aln
    //remove_outliers
    cout << "profile run completed" << endl;
}

void SQLiteProfiler::get_children(string in_id, vector<string> * in_ids, vector<string> * in_keepids){
    Database conn(db);
    string sql = "SELECT ncbi_id FROM taxonomy WHERE parent_ncbi_id = "+in_id;
    sql += " and name_class = 'scientific name';";
    Query query(conn);
    query.get_result(sql);
//	StoreQueryResult R = query.store();
    while(query.fetch_row()){
	string a;
	a = to_string(query.getstr());
	in_ids->push_back(a);
	in_keepids->push_back(a);
    }
}

vector<string> SQLiteProfiler::get_final_children(string id){
    if(id == "1"){
	cout << "special case as id is root" << endl;
	Database conn(db);
	vector<string> allids;
	string sql = "SELECT name,name_class,ncbi_id FROM taxonomy;";
	Query query(conn);
	query.get_result(sql);
	//StoreQueryResult R = query.store();
	while(query.fetch_row()){
	    //string tid = R[j][0].c_str();
	    string tn = query.getstr();
	    string cln = query.getstr();
	    string ncbiid = query.getstr();
	    if(cln.find("scientific")!=string::npos && tn.find("environmental")==string::npos && cln.find("environmental")==string::npos){
		allids.push_back(ncbiid); //was taxon id, now ncbi id
	    }
	}
	query.free_result();
	return allids;
    }else{
	vector<string> ids;
	ids.push_back(id);
	vector<string> keepids;
	//check
	keepids.push_back(id);
	while(ids.size()>0){
	    string tid = ids.back();
	    ids.pop_back();
	    get_children(tid,&ids, &keepids);
	    //cout << ids.size() << endl;
	}
	vector<string> allids;
	vector<string> allnames;

	for(int i=0;i<keepids.size();i++){
	    Database conn(db);
	    string sql = "SELECT name FROM taxonomy WHERE ncbi_id = "+keepids[i];
	    sql += " and name_class = 'scientific name';";
	    Query query(conn);
	    query.get_result(sql);
	    //StoreQueryResult R = query.store();
	    while(query.fetch_row()){
		string name;
		string cl;
		id = keepids[i];
		name = query.getstr();
//			cl = R[j][2].c_str();
//			size_t found1;
//			found1 = cl.find("environmental");
//			size_t found2;
//			found2 = cl.find("scientific");
//			if(found2!=string::npos && found1==string::npos){
		allids.push_back(id);
		allnames.push_back(name);
//			}
	    }
	}
	allids.push_back(id);
	return allids;
    }
}

/*
 * replaces get_final_children and get_children with this
 * that uses the left and right values  to get all the children
 */

vector<string> SQLiteProfiler::get_left_right_children(string id){
    vector<string> allids;
    Database conn(db);
    string sql2 = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = "+id;
    sql2 += " and name_class = 'scientific name';";
    Query query(conn);
    query.get_result(sql2);
    //StoreQueryResult R = query.store();
    string left,right;
    while(query.fetch_row()){
	left = to_string(query.getval());
	right = to_string(query.getval());
    }

    string sql = "SELECT ncbi_id FROM taxonomy WHERE left_value >= "+left;
    sql+="AND right_value <= ";
    sql+=right;
    sql+=" and name_class = 'scientific name' ORDER BY left_value DESC;";
    Query query2(conn);
    query2.get_result(sql);
//	StoreQueryResult R2 = query2.store();
    while(query2.fetch_row()){
	allids.push_back(to_string(query2.getval()));
    }
    return allids;
}

string SQLiteProfiler::get_right_one(vector<string> allids,Query & res){
    while (res.fetch_row()){
	int fd;
	string x;
	x = res.getstr();
	fd = count(allids.begin(),allids.end(),x);
	if (fd > 0)
	    return x;
    }
}

//ncbi one
void SQLiteProfiler::create_distances(string clade_name,map<int, map<int,double> > * numlist){
    numlist->clear();
    //get id for clade name
    // Make SQL string and execute it
    Database conn(db);
//	StoreQueryResult R;
    Query query(conn);
    string sql;
    if (automated == false){
	sql = "SELECT ncbi_id FROM taxonomy WHERE name = '"+clade_name+"' and name_class = 'scientific name';";
	query.get_result(sql);
    }else if (automated == true){
	sql = "SELECT ncbi_id FROM taxonomy WHERE ncbi_id = '"+clade_name+"' and name_class = 'scientific name';";
	query.get_result(sql);
    }
    string cladeid;
    while(query.fetch_row()){
	cladeid = to_string(query.getval());
    }
    cout << cladeid << endl;
    /*
     * get left and right for the cladeid
     */
    string sql2 = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = "+cladeid;
    sql2 += " and name_class = 'scientific name';";
    Query query2(conn);
    query2.get_result(sql2);
//	StoreQueryResult R2 = query2.store();
    string cladeleft,claderight;
    while(query2.fetch_row()){
	cladeleft = to_string(query2.getstr());
	claderight = to_string(query2.getstr());
    }
    /*
     * end left and right
     */
    vector<string> allids = get_final_children(cladeid);
    //vector<string> allids = get_left_right_children(cladeid);
    cout << "clade=" << cladeid << endl;
    //get the route to the clade name
//    for(int i=0;i<names.size();i++){
    map<int,string>::iterator it;
    for (it=profile_id_name_map.begin();it!=profile_id_name_map.end();it++){
	Database conn(db);
	//change to ncbi id
	sql = "SELECT ncbi_id FROM taxonomy WHERE ncbi_id = "+(*it).second+";";//names[i]+";";
//	cout << sql << endl;
	Query query2(conn);
	query2.get_result(sql);
	string nameid = get_right_one(allids, query2);
	sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = "+nameid+" and name_class = 'scientific name';";
	Query query4(conn);
	query4.get_result(sql);
	string parentid = get_right_one(allids,query4);
	vector<string> route;
	while(parentid != cladeid){
	    route.push_back(parentid);
	    nameid = parentid;
	    sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = "+nameid+" and name_class = 'scientific name';";
	    //cout << sql << endl;
	    Query query5(conn);
	    query5.get_result(sql);
	    parentid = get_right_one(allids,query5);
	}
	route.push_back(parentid);

	map<int,double> tdistance;
	map<int,string>::iterator it2;
//	for(int j=0;j<names.size();j++){
	for(it2=profile_id_name_map.begin();it2!=profile_id_name_map.end();it2++){
	    if(it2!=it){
		//sql = "SELECT ncbi_id FROM taxonomy WHERE name = '"+names[j]+"'  and name_class = 'scientific name';";
                //change to ncbi id
		sql = "SELECT ncbi_id FROM taxonomy WHERE ncbi_id = "+(*it2).second+";";//names[j]+";";
		Query query5(conn);
		query5.get_result(sql);
		string jnameid = get_right_one(allids,query5);
		double distance = 0;
		while((int)count(route.begin(),route.end(),jnameid) == 0){
		    distance += 1;
		    sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = "+jnameid+" and name_class = 'scientific name';";
		    Query query6(conn);
		    query6.get_result(sql);
		    jnameid = get_right_one(allids,query6);
		}
		for(int k=0;k<route.size();k++){
		    if(route[k] == jnameid)
			distance += k;
		}
		tdistance[(*it2).first] = distance;
	    }else{
		tdistance[(*it2).first] = SIXES;
	    }
	}
	(*numlist)[(*it).first] = tdistance;
	cout << "distances complete: "<< (*it).second << " " <<(*it).first << endl;
    }
}

//user tree one
void SQLiteProfiler::create_distances_user_tree(vector<string> ifile_names,map<string,string> * numnames
						,map<string,string> * namesnum, map<int, map<int,double> > * numlist){
    //get the list of nodes for which distances are required
    vector<Node *> nodesfordist;
    for(int i=0;i<ifile_names.size();i++){
	for(int j=0;j<userguidetree->getNodeCount();j++){
	    if (userguidetree->getNode(j)->getName()==ifile_names[i])
		nodesfordist.push_back(userguidetree->getNode(j));
	}
    }
    for(int i=0;i<nodesfordist.size();i++){
	vector<double> tdistance;
	for(int j=0;j<nodesfordist.size();j++){
	    if (i == j){
		tdistance.push_back(SIXES);
	    }else{
		double distance = get_distance_between_two_nodes(userguidetree,nodesfordist[i],nodesfordist[j]);
		tdistance.push_back(distance);
	    }
	}
	std::ostringstream stm;
	stm << i;
	numnames->insert( pair<string,string>(stm.str(),nodesfordist[i]->getName()) );
	namesnum->insert( pair<string,string>(nodesfordist[i]->getName(),stm.str()) );
//	numlist->push_back(tdistance);
    }
}

void SQLiteProfiler::get_shortest_distance_with_dicts(
                vector<int> & nums,
				map<int, map<int, double> > & numlist,
                int * sd_aln1,
                vector<int> * sd_alns2) {

    /* this function finds the alignment pairs separated by the shortest taxonomic distance  */
    
    map<int,double> distances;
    double shortestdistance = 10000;

    for(int i = 0; i < nums.size(); i++) {
        bool keepD = false;
        for(int j=0; j < first_profiles_dbids.size(); j++) {
            double distance = numlist[nums[i]][first_profiles_dbids[j]];
                if(distance < shortestdistance && distance != SIXES) {
                shortestdistance = distance;
                *sd_aln1 = nums[i];
                keepD = true;
            }
        }
        if (keepD == true) {
            distances = numlist[nums[i]];
        }
    }

    sd_alns2->clear();
    for(int j = 0; j < first_profiles_dbids.size(); j++) {
        if(distances[first_profiles_dbids[j]] == shortestdistance && distances[first_profiles_dbids[j]] != SIXES) {
            sd_alns2->push_back(first_profiles_dbids[j]); // should this be nums
        //	    cout << "f " <<first_profiles_dbids[j] << " " <<distances[first_profiles_dbids[j]] << endl;
        }
    }
}

void SQLiteProfiler::clean_before_profile(string filename) {

    /* just cleaning the file that is there
     * 0.1 is the current limit
     * percent should be, if you are missing more than this, then remove, so 0.9 is now and 0.5 is genome removal
     */

    double threshold = 0.9;

    // read the alignment into tempalseqs
    FastaUtil fu;
    vector<Sequence> tempalseqs;
    fu.readFile(profilefoldername + filename, tempalseqs);

    cout << "cleaning seqs for " << filename << endl;
    int seqlength = tempalseqs[0].get_sequence().size();
    float nseqs = float(tempalseqs.size());
    vector<int> cols_to_remove;

    // for each column in the alignment
    for(int j = 0; j < seqlength; j++) {
        int gaps = 0;
        char cell;
        // get the number of sequences with missing data
        for(int i = 0; i < tempalseqs.size(); i++) {
            cell = tempalseqs[i].get_sequence()[j];
            if(cell == '-' || cell == 'N' || cell == 'n')
                gaps += 1;
        }
        // if we exceed the threshold, mark this column for removal
        double prop_missing = gaps/nseqs;
        if (prop_missing > threshold)
            cols_to_remove.push_back(j);
    }

    // for each sequence
    for(int i = 0; i < tempalseqs.size(); i++) {
        string a;
        // for each column
        for (int j = 0; j < seqlength; j++) {
            // record all columns not in the exclusion list
            if (count(cols_to_remove.begin(), cols_to_remove.end(), j) == 0)
                a += tempalseqs[i].get_sequence()[j];
        }
        // record cleaned sequence
        tempalseqs[i].set_sequence(a);
    }

    // remove the input (dirty) file, replace it with the clean one
    remove((profilefoldername + filename).c_str());
    fu.writeFileFromVector(profilefoldername + filename, tempalseqs);
}

void SQLiteProfiler::clean_dbseqs(int alignid){
    double percent = 0.9;
    vector<Sequence> tempalseqs;
    gene_db.get_profile_align_seqs(alignid,tempalseqs);
    cout << "cleaning seqs" << endl;
    int seqlength = tempalseqs[0].get_aligned_seq().size();
    float fseql = float(tempalseqs.size());
    vector<int> removeem;
    for(int j=0;j<seqlength;j++){
	int gaps = 0;
	for(int i=0;i<tempalseqs.size();i++){
	    if(tempalseqs[i].get_aligned_seq()[j] == '-' || tempalseqs[i].get_aligned_seq()[j] == 'N' || tempalseqs[i].get_aligned_seq()[j] == 'n')
		gaps += 1;
	}
	double curp = gaps/fseql;
	if (curp > percent){
	    removeem.push_back(j);
	}
    }
    for(int i=0;i<tempalseqs.size();i++){
	string a;
	for (int j=0;j<seqlength;j++){
	    if (count(removeem.begin(),removeem.end(),j)==0)
		a += tempalseqs[i].get_aligned_seq()[j];
	}
	tempalseqs[i].set_aligned_seq(a);
    }
    gene_db.update_profile_align_seqs(alignid,tempalseqs);
}

/*
 * TODO: fix the orphans
 */
int SQLiteProfiler::profile(map<int, map<int,double> > numlist) {

    /* this processes through numlist (a list of database alignment ids), profile-aligning them to
     * one another until all alignments have been profiled.
     * 
     * we keep track of original and profile alignments in two vectors: original_alns_to_profile and
     * profile_alns_to_profile, which store the alignments' database ids. we process through these,
     * deleting the alignment ids (from the respective vector) as they are combined into a new profile
     * alignment, and adding the id of the new profile alignment to the profile_alns_to_profile
     * list. we stop when original_alns_to_profile is empty and only one profile id remains; the
     * corresponding final alignment contains all the original alignments referenced in numlist. */

    cout << "writing everything to record.log" << endl; // really? this may not actually be the case...
    string recordname = profilefoldername + "/record.log"; 
    ofstream ofs(recordname.c_str());

    bool muscle = true;

    // initiate arrays to store alignment ids to profile.
    // all ids are database row ids from the profile_alignments table.
    vector<int> profile_alns_to_profile;

    // the original alignments have been copied to the profile_alignments table in the database.
    // this array stores their row ids in that table.
    vector<int> original_alns_to_profile(first_profiles_dbids.begin(), first_profiles_dbids.end());

    // DEBUG: logging output
//    cout << "dbids for starting profiles (original alignments): " << endl;
//    for(int i = 0; i < original_alns_to_profile.size(); i++) {
//        cout << original_alns_to_profile[i] << endl;}

//    exit(0);

    if (updatedb == true) {

        // reduce the alignments to profile to just those that have been updated
        original_alns_to_profile = updatednums;
        cout << "original alignments remaining to profile: ";
        for (int i = 0; i < original_alns_to_profile.size(); i++)
            cout << original_alns_to_profile[i] << " ";
        cout << endl;

        // add the existing profiles to profile_alns_to_profile
        profile_alns_to_profile = updatedprofiles;
        cout << "intermediate profile alignments remaining to cross-align: ";
        for (int i = 0; i < profile_alns_to_profile.size(); i++)
            cout << profile_alns_to_profile[i] << " ";
        cout << endl;
    }

    // use this to flag when to change the output name to FINAL.aln
    int last_profile_file;

    // first, profile all the original alignments
    while (original_alns_to_profile.size() > 0) {

        // find the alignment pairs separated by the shortest taxonomic distance
        // more than one alignment may be equally close to shortest_dist_aln1
        int shortest_dist_aln1;
        vector<int> * shortest_dist_alns2 = new vector<int>();
        get_shortest_distance_with_dicts(
                    original_alns_to_profile,
                    numlist,
                    &shortest_dist_aln1,
                    shortest_dist_alns2);

        // DEBUG: log closest alignment pairs
//        cout << "first alignment to profile: " << shortest_dist_aln1 << endl;
//        cout << "shortest_dist_alns2 = ";
//        for (int i = 0; i < shortest_dist_alns2->size(); i++)
//            cout << shortest_dist_alns2->at(i) << " ";
//        cout << endl;

        // remove the starting alignment from the set of alignments to profile
        vector<int>::iterator it;
        it = find(original_alns_to_profile.begin(), original_alns_to_profile.end(), shortest_dist_aln1);
        original_alns_to_profile.erase(it);

        if(shortest_dist_alns2->size() == 1) {
            // if there is only one closest match to shortest_dist_aln1

            int firstfile = shortest_dist_aln1;
            int secondfile;
            bool sd_aln2_already_profiled = false;

            // see if this match is already within a profile alignment
            secondfile = gene_db.get_deepest_profile_for_alignment(shortest_dist_alns2->at(0));
            if (secondfile != -1)
                // if it is, use the containing profile alignment (secondfile)
                sd_aln2_already_profiled = true;

            // if not, use the original alignment
            if(sd_aln2_already_profiled == false)
                secondfile = shortest_dist_alns2->at(0);

            // create a new profile alignment record, add it to the list of profile alns to cross-align
            int profileout_id = gene_db.add_profile_alignment(firstfile, secondfile);
            profile_alns_to_profile.push_back(profileout_id);

            // do the profile alignment
            int confirmed_profile_id = make_muscle_profile(firstfile, secondfile, profileout_id);
            last_profile_file = confirmed_profile_id;

            // store the new profile alignment in the db
            clean_before_profile("TEMPOUT.PROFILE");
            match_and_add_profile_alignment_to_db(profileout_id);

            // now remove the second input alignment, since it has been stored in the new profile aln
            if(sd_aln2_already_profiled == false) {

                // if input2 was an original alignment, remove it from that array
                vector<int>::iterator it;
                it = find(original_alns_to_profile.begin(), original_alns_to_profile.end(), secondfile);
                original_alns_to_profile.erase(it);

            } else {

                // otherwise remove it from the profiles to cross-profile array
                vector<int>::iterator it;
                it = find(profile_alns_to_profile.begin(), profile_alns_to_profile.end(), secondfile);
                profile_alns_to_profile.erase(it);
            }
        
        } else { // there is more than one equally closest alignment to the starting alignment

            int firstfile = shortest_dist_aln1;
            int secondfile;
            bool sd_aln2_already_profiled = false;

            // walk the alignments in the list of equally closest matches...
            for(int j = 0; j < shortest_dist_alns2->size(); j++) {
                // ...and check if each is a profile alignment
                secondfile = gene_db.get_deepest_profile_for_alignment(shortest_dist_alns2->at(j));

                // DEBUG
//                cout << "Attempting to compare secondfile = " << secondfile << " == -1. Result is: ";

                if (secondfile != -1) {
//                    cout << "FALSE" << endl; // DEBUG
//                    cout << "Found a matching profile alignment with id " << secondfile;
//                    cout << " for sd_aln id " << shortest_dist_alns2->at(j) << endl;
                    // when we find a profile alignment, break
                    sd_aln2_already_profiled = true;
                    break;

                }
                // DEBUG
//                else {
//                    cout << "TRUE" << endl;
//                }
            }

            if (sd_aln2_already_profiled == false) {
                int bestsn;
                double bestscore = 0;
                for(int i = 0; i < shortest_dist_alns2->size(); i++) {
                    make_muscle_profile(firstfile,shortest_dist_alns2->at(i),0); // needs to be a fake one (sas) - huh? (ceh)
                    double score  = get_muscle_spscore("TEMPOUT.PROFILE");

                    if(score > bestscore) {
                        bestscore = score;
                        bestsn = shortest_dist_alns2->at(i);
                    }
            //		cout << "score: " << score << endl;
                }
                  
                secondfile = bestsn;
//                cout << "secondfile: " << secondfile << endl;
                
                int profileout_id = gene_db.add_profile_alignment(firstfile,secondfile);

                // profile
                int outfile_id = make_muscle_profile(firstfile,secondfile,profileout_id);
                clean_before_profile("TEMPOUT.PROFILE");
                match_and_add_profile_alignment_to_db(profileout_id);

                //clean
                profile_alns_to_profile.push_back(profileout_id);
                last_profile_file = outfile_id;

                // DEBUG
//                cout << "looking for secondfile = " << secondfile << " in original_alns_to_profile" << endl;
//                cout << "original_alns_to_profile.begin = " << &(*original_alns_to_profile.begin()) << endl;
//                cout << "original_alns_to_profile.end = " <<  &(*original_alns_to_profile.end()) << endl;

                vector<int>::iterator it;                
                it = find(original_alns_to_profile.begin(), original_alns_to_profile.end(), secondfile);

                // DEBUG
                // dumping data to see if we can identify what is going on with the profile delete mismatch bug
//                cout << "last_profile_id = " << last_profile_file << endl;
//                cout << "profile_alns_to_profile[profile_alns_to_profile.size() - 1] = " << profile_alns_to_profile[profile_alns_to_profile.size() - 1] << endl;
//                cout << "original_alns_to_profile: " << endl;
//                for(int i = 0; i < original_alns_to_profile.size(); i++) {
//                    cout << original_alns_to_profile[i] << endl;} 
//                cout << "profile_alns_to_profile: " << endl;
//                for(int i = 0; i < profile_alns_to_profile.size(); i++) {
//                    cout << profile_alns_to_profile[i] << endl;}
//                cout << "Attempting to delete id " << *it << " from original_alns_to_profile" << endl;

                // remove the child alignment from from the table; it is in a profile now
                original_alns_to_profile.erase(it);
            
            } else {
            int profileout_id = gene_db.add_profile_alignment(firstfile,secondfile);
            int s2 = make_muscle_profile(firstfile,secondfile,profileout_id);
            clean_before_profile("TEMPOUT.PROFILE");
            match_and_add_profile_alignment_to_db(profileout_id);
            profile_alns_to_profile.push_back(profileout_id);
            last_profile_file = s2;
            vector<int>::iterator it;
            it = find (profile_alns_to_profile.begin(), profile_alns_to_profile.end(), secondfile);
            profile_alns_to_profile.erase(it);
            }
        }
        delete shortest_dist_alns2;
    //	exit(0);
    }
    cout << "processing profile files" << endl;
    for(int i=0; i < profile_alns_to_profile.size(); i++){
        cout << profile_alns_to_profile[i] << endl;}


    while(profile_alns_to_profile.size() > 1){
//TODO: choose the best, for now just choose any two
	int firstfile = profile_alns_to_profile[0];
	int secondfile = profile_alns_to_profile[1];
	vector<int>::iterator it;
	it = find (profile_alns_to_profile.begin(), profile_alns_to_profile.end(), firstfile);
	profile_alns_to_profile.erase(it);
	it = find (profile_alns_to_profile.begin(), profile_alns_to_profile.end(), secondfile);
	profile_alns_to_profile.erase(it);
	
	int profileout_id = gene_db.add_profile_alignment(firstfile,secondfile);
	int s2 = make_muscle_profile(firstfile,secondfile,profileout_id);
	clean_before_profile("TEMPOUT.PROFILE");
	match_and_add_profile_alignment_to_db(profileout_id);
	profile_alns_to_profile.push_back(s2);
	last_profile_file = s2;
    }

    return last_profile_file;
}

int SQLiteProfiler::make_muscle_profile(int profile1,int profile2,int outfile_id){

    /* accepts the database ids for two preexisiting profile alignments, writes
     * these alignments to files, and calls the external program muscle to align
     * them together (profile alignment). the resulting alignment is stored in a
     * file named TEMPOUT.PROFILE. if this file is created successfully, then the
     * third function argument, representing the next database id (the one to be
     * used for the new profile alignment) is returned. */
    
    // delete any previous temp alignment files
    remove((profilefoldername+"TEMP1.PROFILE").c_str());
    remove((profilefoldername+"TEMP2.PROFILE").c_str());
    remove((profilefoldername+"TEMPOUT.PROFILE").c_str());

    // write the alignments to temp files
    gene_db.write_profile_alignment_to_file(profile1, profilefoldername+"TEMP1.profile");
    gene_db.write_profile_alignment_to_file(profile2, profilefoldername+"TEMP2.profile");

    // build the muscle command call
    string cmd = "muscle -profile -maxmb 5000 -in1 ";
    cmd += profilefoldername + "TEMP1.profile -in2 ";
    cmd += profilefoldername + "TEMP2.profile -out ";
    cmd += profilefoldername + "TEMPOUT.PROFILE 2> ";
	cmd += profilefoldername + "muscle.out";

    cout << "aligning" << endl;
    cout << cmd << endl;

    // call muscle to do profile alignment
    system(cmd.c_str());

    // if new profile alignment was not created, test will fail and program will exit
    test_outfile_exists(profilefoldername+"TEMPOUT.PROFILE");

    // otherwise return the new profile's id, passed in when make_muscle_profile was called
    return outfile_id;
}


string SQLiteProfiler::get_name_from_tax_id(string taxid){
    Database conn(db);
    string sql = "SELECT name FROM taxonomy WHERE ncbi_id = "+taxid+" AND name_class = 'scientific name'";
    Query query(conn);
    query.get_result(sql);
//	StoreQueryResult R = query.store();
    string ret;
    while(query.fetch_row()){
	string a;
	a = query.getstr();
	ret =a;
    }
    return ret;
}

/*
 * rename the FINAL.aln.cln file to FINAL.aln.cln.rn using the table in ITS.gi
 */
void SQLiteProfiler::rename_final_alignment(int alignid){
    string fn1 = gene_name+".FINAL.aln.rn";
    gene_db.write_profile_alignment_with_names_to_file(alignid,fn1,false);
    fn1 = gene_name+".FINAL.aln";
    gene_db.write_profile_alignment_with_names_to_file(alignid,fn1,true);
}

double SQLiteProfiler::get_muscle_spscore(string filename){
    remove("prolog");
    string cmd = "muscle -maxmb 5000 -spscore ";
    cmd += profilefoldername;
    cmd += filename;
    cmd += " -log ";
	cmd += profilefoldername+"prlog 2> ";
	cmd += profilefoldername+"muscle.out";
//    cout << cmd << endl;
    FILE *fp = popen(cmd.c_str(), "r" );
    char buff[1000];
    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	string line(buff);
    }
    pclose( fp );
    //read file
    double score = 0;
    ifstream ifs((profilefoldername+"prlog").c_str());
    string line;
    while(getline(ifs,line)){
	TrimSpaces(line);
	size_t se = line.find("="); // Find the first character position from reverse af
	// if all spaces or empty return an empty string
	if (se != string::npos){
	    vector<string> tokens;
	    string del("=");
	    Tokenize(line, tokens, del);
	    score = atof(tokens[2].c_str());
	}
    }
    ifs.close();
    return score;
}

void SQLiteProfiler::test_outfile_exists(string filename){
    struct stat filestatus;
    stat(filename.c_str(), &filestatus );
    cout << filestatus.st_size << " bytes\n";
    if (int(filestatus.st_size) <= 64){
	cerr << "problem: empty file "<< filename<< " created" << endl;
	exit(0);
    }
}

/*
 * assumes that the id of the seq is the sqlite id 
 */
void SQLiteProfiler::match_and_add_profile_alignment_to_db(int profileid){
    FastaUtil fu;
    vector<Sequence> tempalseqs;
    fu.readFile(profilefoldername+"TEMPOUT.PROFILE",tempalseqs);
    vector<Sequence> seqs;
    for(int i=0;i<tempalseqs.size();i++){
	Sequence tseq(tempalseqs[i].get_id(),tempalseqs[i].get_sequence());
	seqs.push_back(tseq);
    }
    cout << "adding sequences to profile alignment table" << endl;
    gene_db.add_sequences_for_profile_alignment(profileid,seqs);
}



//THIS WAS THE OLD quicktree based outlier removal
/*
void SQLiteProfiler::calculate_for_removal_quicktree(vector<Sequence> * seqs,
		map<string,double> & allmeans){
    FastaUtil seqwriter;
    const string fn1 = "TEMPFILES/tempremoval";
    seqwriter.writeFileFromVector(fn1,*seqs);

    string phcmd = "phyutility -concat -in TEMPFILES/tempremoval -out TEMPFILES/outfile.nex";
    FILE *phfp = popen(phcmd.c_str(), "r" );
    pclose( phfp );

    cout << phcmd << endl;

    ifstream infile;
    ofstream outfile;
    infile.open ("TEMPFILES/outfile.nex",ios::in);
    outfile.open ("TEMPFILES/outfile.stoc",ios::out);
    bool begin = false;
    bool end = false;
    string line;
    // convert to stockholm format
    while(getline(infile,line)){
	if (line.find("MATRIX") != string::npos){
	    begin = true;
	}else if ((begin == true && end == false) && line.find_first_of(";") != string::npos){
	    end = true;
	}else if (begin == true && end == false){
	    std::string::size_type begin = line.find_first_not_of("\t");
	    //std::string::size_type end   = line.find_last_not_of("\t");
	    std::string::size_type end = line.size();
	    std::string trimmed = line.substr(begin, end-begin + 1);
	    outfile << trimmed << endl;
	}
    }
    infile.close();
    outfile.close();

    const char * cmd = "quicktree -in a -out m TEMPFILES/outfile.stoc > TEMPFILES/dist";
    cout << "calculating distance" << endl;
    FILE *fp = popen(cmd, "r" );
    char buff[1000];
    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	string line(buff);
    }
    pclose( fp );

    //new
    vector<string> ids;
    vector<vector<double> > nums;
    bool first = true;

    ifstream pfile ("TEMPFILES/dist");
    vector<string> tokens;
    int nspecies;
    if (pfile.is_open()){
	int curspecies = 0;
	while (! pfile.eof() ){
	    // record the id as the first token
	    // and the numbers below as the numbers
	    if (first == false){
		getline (pfile,line);
		string del(" \t");
		tokens.clear();
		Tokenize(line, tokens, del);
		if(tokens.size() >= 1){
		    ids[curspecies] = tokens.at(0);
		    double n1;
		    for(int j = curspecies; j < nspecies;j++){
			n1 = atof(tokens.at(j+1).c_str());
			nums[curspecies][j] = n1;
			nums[j][curspecies] = n1;
		    }
		    curspecies += 1;
		}
	    }else{
		first = false;
		getline (pfile,line);
		TrimSpaces(line);
		nspecies = atoi(line.c_str());
		vector<double> cols(nspecies, 0);
		nums = vector< vector<double> >(nspecies, cols);
		ids = vector<string>(nspecies);
	    }
	}
	pfile.close();
    }
// calculate the means
    vector<double> mns(ids.size());
    for(int i=0;i<nums.size();i++){
	allmeans[ids[i]] = mean(nums[i]);
    }
}

void SQLiteProfiler::remove_outliers(){
    FastaUtil seqreader;
    vector<Sequence> * sequences = new vector<Sequence>();
    seqreader.readFile(profilefoldername+"FINAL.aln", *sequences);
    int numseqs = sequences->size();
    cout << numseqs << endl;
    map<string,double> allmeans;
    if(numseqs < 5000){
	calculate_for_removal_quicktree(sequences,allmeans);
    }else{
	int NBREAKS = 10;
	for (int i=0;i<NBREAKS;i++){
	    vector<Sequence> tempsc1;
	    if((i+1) < NBREAKS){
		for(int j=(i*(numseqs/NBREAKS));j < ((numseqs/NBREAKS)*(i+1));j++){
		    //TODO : there was a pointer problem here
		    tempsc1.push_back(sequences->at(j));
		}
	    }else{
		for(int j=(i*(numseqs/NBREAKS));j < numseqs;j++){
		    //TODO : there was a pointer problem here
		    tempsc1.push_back(sequences->at(j));
		}
	    }
	    calculate_for_removal_quicktree(&tempsc1,allmeans);
	}
    }
    // calculate the means
    vector<double> mns(numseqs);
    for(int i=0;i<numseqs;i++){
	//TODO : there was a pointer problem here
	mns[i] = allmeans[sequences->at(i).get_id()];
    }
    double sd = stdev(mns);
    double mn = mean(mns);
    double dev = sd*1.5+mn;//obviously changeable
    // read in the fasta file
    vector<Sequence> sc1;
    for(int i=0;i<mns.size();i++){
	if (mns[i] < dev){
	    sc1.push_back(sequences->at(i));
	}else{
	    //cout << i << endl;
	}
    }
    FastaUtil seqwriter;
    const string fn1 = profilefoldername+"FINAL.aln.cln";
    seqwriter.writeFileFromVector(fn1,sc1);
}
*/

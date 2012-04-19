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
 * SQLiteConstructor.cpp
 */



#include "SWPS3_matrix.h"

#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <limits>
#include <set>

using namespace std;

#include "omp.h"
#include "libsqlitewrapped.h"
#include "tree.h"
#include "tree_reader.h"
#include "sequence.h"
#include "fasta_util.h"
#include "SQLiteConstructor.h"
#include "utils.h"

//the number for the mafft alignment sample
#define RANDNUM 1000
#define ALIGNLIMIT 7000

//public

template <class T>
inline std::string to_string (const T& t){
    std::stringstream ss;
    ss << t;
    return ss.str();
}

SQLiteConstructor::SQLiteConstructor(string cn, vector <string> searchstr, string genen, string genedb,
             double mad_cut,double cover, double ident, string dbs, 
             string known_seq_filen, bool its, int numt,bool autom,
             bool inupdatedb, string inupdatefile): clade_name(cn),search(searchstr),
                      gene_name(genen), gene_db_name(genedb),
                      mad_cutoff(mad_cut),coverage(cover),
                      identity(ident),db(dbs),useITS(its),
                      numthreads(numt),automated(autom),
                      updateDB(inupdatedb){
    FastaUtil seqReader;
    //added updating seqs from db
    //added updating seqs from file
    if(inupdatefile.length() > 0){
	updateFILE = true;
	updatef = inupdatefile;
    }else{
	updateFILE = false;
    }
    //initialize gene database
    gene_db= GeneDB(gene_db_name);
    genefoldername = genen+"_TEMPFILES/";
    known_seqs = new vector<Sequence>();
    user_seqs = new vector<Sequence>();
    seqReader.readFile(known_seq_filen, *known_seqs);
    ncbi_saturation = true;
    userskipsearchdb = false;
    skipdbcheck = false;
    justseqquery = false;
}

/*
 * This will supercede the previous one in the main file. If you have higher taxa
 * in the listfile and you say they are wildcards then it will supercede if you
 * have higher taxa and don't want to pick one as above.
 */
void SQLiteConstructor::set_only_names_from_file(string filename, bool containshi, bool containswi){
    onlynamesfromfile = true;
    containshigher = containshi;
    containswild = containswi;
    if (containswild == true)
	containshigher = false;
    listfilename = filename;
}

void SQLiteConstructor::set_exclude_names_from_file(string filename){
    excludenamesfromfile = true;
    excludefilename = filename;
}

void SQLiteConstructor::set_exclude_gi_from_file(string filename){
    excludegifromfile = true;
    exclude_gi_filename = filename;
}

void SQLiteConstructor::set_include_gi_from_file(string filename){
    includegifromfile = true;
    include_gi_filename = filename;
}

void SQLiteConstructor::set_user_skip_search(){
    userskipsearchdb = true;
}

void SQLiteConstructor::set_justseqquery(bool setit){
    justseqquery = setit;
}

/*
 * if the user guide tree covers less than 50% of the taxa, falling back
 * on NCBI
 * this will allow for skipping the check to the ncbi db, not really that 
 * useful, but good for things like simulated data
 */
void SQLiteConstructor::set_user_guide_tree(string filename,bool inskipdbcheck){
    usertree = true;
    usertreefile = filename;
    //read the tree here
    TreeReader nw;
    ifstream infile2(filename.c_str());
    if (!infile2){
	cerr <<"user guide tree: "<<filename<< " doesn't exist"<<endl;
	exit(0);
    }
    vector<string> lines;
    string line;
    while(getline(infile2,line)){
	lines.push_back(line);
    }
    infile2.close();
    userguidetree = nw.readTree(lines[0]);
    cout << "read user guide tree with : "<<userguidetree->getExternalNodeCount() << " tips"<<endl;
    //need to number the internal nodes so that they can be used for the run 
    int count = 0;
    for(int i=0;i<userguidetree->getInternalNodeCount();i++){
	userguidetree->getInternalNode(i)->setName(to_string(count));
	count +=1;
    }
    //procedure would be to 
    //get the ncbi taxa for each of these
    //associating that information with each of the tip nodes as comment
    skipdbcheck = inskipdbcheck;
    if (skipdbcheck==false){
	Database conn(db);
	//TODO: need to make this faster
	cout << "matching user guide tree names: "<<userguidetree->getExternalNodeCount() << endl;
	for (int i=0;i<userguidetree->getExternalNodeCount();i++){
	    string tname = userguidetree->getExternalNode(i)->getName();
	    //search for the name in the database
	    string sql = "SELECT ncbi_id FROM taxonomy WHERE edited_name = '"+tname+"' OR ncbi_id = '"+tname+"';";
	    Query query(conn);
	    query.get_result(sql);
	    int count1 = 0;
	    bool nameset = false;
	    string nameval;
	    while(query.fetch_row()){
		nameset = true;
		nameval = to_string(query.getval());
		count1+=1;
	    }
	    query.free_result();
	    if (nameset==true){
		userguidetree->getExternalNode(i)->setComment(nameval);
		//cout << tname << "="<<nameval<<endl;
	    }else{
		cerr <<tname << " is not in the ncbi database as a number or name"<<endl; 
	    }
	}
    }
    //this will change the search to use the user tree instead of the ncbi tree
    ncbi_saturation = false;
}

/*
 * allows for adding user seqs
 * will also allow for skipping a db check. really not that useful unless
 * there is simulated data or something similar
 */
void SQLiteConstructor::set_user_fasta_file(string filename,bool skipdbcheck){
    userfasta = true;
    userfastafile = filename;
    FastaUtil fu;
    cout <<"reading user fasta file: "<<userfastafile <<endl;
    fu.readFile(userfastafile, *user_seqs);
    cout << "successfully read " << user_seqs->size()<< " user seqs"<<endl;
    //need to edit the names if there aren't already changed
    //should go into the log
    for (int i=0;i<user_seqs->size();i++){
	string tstring (user_seqs->at(i).get_id());
	fix_bad_chars_for_seq_names(tstring);
	//if user is not in the front then add it
	if (tstring.find("user_") != 0)
	  tstring = "user_"+tstring;
	cout <<"changing "<< user_seqs->at(i).get_id()<<" -> "<< tstring << endl;
	user_seqs->at(i).set_id(tstring);
    }
    //if the ids can be found in the database, this will go ahead and make that link
    //this should be able to link to the genedb eventually as well
    //store the ncbi id in the comments
    if(skipdbcheck==false){
	Database conn(db);
	//TODO: need to make this faster
	cout << "trying to match sequence names to ncbi name"<< endl;
	for (int i=0;i<user_seqs->size();i++){
	  string tname = user_seqs->at(i).get_id().substr(5,user_seqs->at(i).get_id().size());
	    
	    //search for the name in the database
	    string sql = "SELECT ncbi_id FROM taxonomy WHERE edited_name = '"+tname+"' OR ncbi_id = '"+tname+"';";
	    Query query(conn);
	    query.get_result(sql);
	    int count1 = 0;
	    bool nameset = false;
	    string nameval;
	    while(query.fetch_row()){
		nameset = true;
		nameval = to_string(query.getval());
		count1+=1;
	    }
	    query.free_result();
	    if (nameset==true){
		user_seqs->at(i).set_ncbi_tax_id(nameval);
		cout << tname << "="<<nameval<<endl;
	    }else{
		cerr <<tname << " is not in the ncbi database as a number or name"<<endl;
		user_seqs->at(i).set_ncbi_tax_id("0");
	    }
	}
    }
    if(usertree == true){
	cout << "matching user fasta seqs to user guide tree" << endl;
	int count = 0;
	for(int i=0;i<userguidetree->getExternalNodeCount();i++){
	    string tname = userguidetree->getExternalNode(i)->getName();
	    //cout << tname << endl;
	    for(int j=0;j<user_seqs->size();j++){
	      if (tname == user_seqs->at(j).get_id() || tname == user_seqs->at(j).get_comment() || ("user_"+tname == user_seqs->at(j).get_id()) || ("user_"+tname == user_seqs->at(j).get_comment())){
		    user_fasta_node_map[&(user_seqs->at(j))] = userguidetree->getExternalNode(i);
		    count +=1;
		}
	    }
	}
	cout << "matches: "<<count<< " prop:"<< count/user_seqs->size() << endl;
    }
    cout <<"finished reading user fasta file" << endl;
}

int SQLiteConstructor::run(){
    string logn = gene_name+".log";
    logfile.open(logn.c_str());
    string gin = gene_name+".gi";
    string uffa = gene_name+".userfasta";
    //updatedb code
    vector<Sequence> stored_seqs;
    vector<string> stored_seqs_ncbis;
    vector<string> stored_user_seqnames;
    write_EDNAFILE();//if it doesn't exist
    if (updateDB == true)
	gene_db.initialize(false);
    else
	gene_db.initialize(true);
    if(updateDB == true){
	cout << "processing existing seqs" << endl;
	gene_db.get_all_sequences(stored_seqs);
	for (int i=0;i<stored_seqs.size();i++){
	    if(stored_seqs[i].get_name().find("user_") != 0)//not a user seq with no id
		stored_seqs_ncbis.push_back(stored_seqs[i].get_id());
	    else//user seq
		stored_user_seqnames.push_back(stored_seqs[i].get_id());
	}
	cout << "existing seqs: " << stored_seqs.size() << endl;
	cout << "existing db: " << stored_seqs_ncbis.size() << " user: " << stored_user_seqnames.size() << endl; 
    }else{
	gifile.open(gin.c_str(),fstream::out);
	gifile << "ncbi_tax_id\tgi_number\tedited_name" << endl;
	//output the user fasta seqs
	if(userfasta){
            ufafile.open(uffa.c_str(),fstream::out);
	    ufafile << "edited_name\tncbi_taxon_id"<<endl;
	}
    }

    mkdir(genefoldername.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH);

    vector<vector<string> > start_res;
    vector<Sequence> startseqs; 
    Database conn(db);
    string sname_id;
    if(userskipsearchdb == false){
	first_seq_search_for_gene_left_right(start_res);

	//make connection to database
	cout << "connected to " << db << endl;
	vector<int> R;
	//just some edits to use the ncbi ids instead of the names for automated runs
	if (automated == false){
	    Query query(conn);
	    query.get_result("SELECT ncbi_id FROM taxonomy WHERE name = '"+clade_name+"'");
	    while(query.fetch_row()){
		R.push_back(query.getval());
	    }
	    query.free_result();
	}else if(automated == true){
	    Query query(conn);
	    query.get_result("SELECT ncbi_id FROM taxonomy WHERE ncbi_id = "+clade_name);
	    while(query.fetch_row()){
		R.push_back(query.getval());
	    }
	    query.free_result();
	}
	//start with a name -- get the broad name clade id
	cout << "Found " << R.size()<< " taxon ids:" << endl;
	for(int i = 0; i < R.size(); ++i){
	    cout << R[i] << endl;
	}
	int name_id;
	//string sname_id;
	name_id = R[0];
	sname_id = to_string(R[0]).c_str();
	cout << "Will be using " << name_id << endl;

	//start with a set of seqs given the first clade name and the regions
	startseqs = first_get_seqs_for_name_use_left_right(name_id, start_res);

	cout << "first: " << startseqs.size() << endl;

	//if excluding gi's from file
	if(excludegifromfile == true){
	    startseqs = exclude_gis_from_file(startseqs);
	    cout << "after excluding gis: " << startseqs.size() << endl;
	}

	//if including gi's from file
	if(includegifromfile == true){
	    startseqs = include_gis_from_file(startseqs);
	    cout << "after including gis: " << startseqs.size() << endl;
	}

	//if only using the names from a list
	if(onlynamesfromfile == true){
	    startseqs = use_only_names_from_file(startseqs);
	    cout << "after names: " << startseqs.size() << endl;
	}

	//if excluding names from file

	if(excludenamesfromfile == true){
	    startseqs = exclude_names_from_file(startseqs);
	    cout << "after excluding names: " << startseqs.size() << endl;
	}

	if (startseqs.size() == 0){
	    cout << "there were no seqs left after blast" << endl;
	    exit(0);
	}
    }//end skipseach==false

    //use blast to idenify seqs and rc
    vector<Sequence> * keep_seqs = new vector<Sequence>();
    if(justseqquery == true)
	get_same_seqs_openmp_SWPS3_justquery(startseqs,keep_seqs);
    else 
	get_same_seqs_openmp_SWPS3(startseqs,keep_seqs);
    cout << "blasted: "<< keep_seqs->size() << endl;
    if (justseqquery == true){
	cout << "scores written to "<< gene_name<<".seqquery" << endl;
	exit(0);
    }
    //assuming for now that all the user seqs are hits

    remove_duplicates_SWPS3(keep_seqs);
    cout << "dups: "<< keep_seqs->size() << endl;
    //if userguidetree overlaps with less than a certain percentage, usertree = false
    if(usertree == true && userskipsearchdb == false){
	double overlap = get_usertree_keepseq_overlap(keep_seqs);
	overlap = 0.5;//TODO: take this out
	if (overlap < 0.5){
	    usertree = false;
	    userguidetree = NULL;
	    cerr << "user guide has only "<< overlap<< " overlap so exiting"<<endl;
	    cout << "comment userguidetree with # in the config file to skip"<<endl;
	    exit(0);
	}
    }

    reduce_genomes(keep_seqs);
	
    //get the list of files and record the higher taxa name and 
    //add the additional sequences to the right hierarchy
    vector<string> sname_ids;
    vector<string> snames;
    if(updateDB == true){
	vector<int> toremove;
	for(int j=0;j<keep_seqs->size();j++){
	    if(count(stored_seqs_ncbis.begin(),stored_seqs_ncbis.end(),keep_seqs->at(j).get_ncbi_tax_id()) > 0){
		toremove.push_back(j);
	    }
	}
	for(int j=0;j<toremove.size();j++){
	    keep_seqs->erase(keep_seqs->begin()+toremove[toremove.size()-(j+1)]);
	}
	//end the program if there is nothing to update
	if(keep_seqs->size() == 0 ){
	    cout << "There are no new DB sequences to add." << endl;
	    gifile.close();
	}
	cout << "total size of seqs to add (update):" << keep_seqs->size() << endl;
	for(int j=0;j<keep_seqs->size();j++){
	    cout << "adding: " << keep_seqs->at(j).get_ncbi_tax_id() << endl;
	}
	gene_db.add_seqs_to_db(keep_seqs);
	//add the right gi numbers before add the rest of the seqs are added to keep_seqs
	write_gi_numbers(keep_seqs);
	gifile.close();
	//adding the update for user seqs
	toremove.clear();
	for(int j=0;j<user_seqs->size();j++){
	  if(count(stored_user_seqnames.begin(),stored_user_seqnames.end(),user_seqs->at(j).get_id()) > 0){
		toremove.push_back(j);
	    }
	}
	for(int j=0;j<toremove.size();j++){
	    user_seqs->erase(user_seqs->begin()+toremove[toremove.size()-(j+1)]);
	}
	if(user_seqs->size() == 0 ){
	    cout << "There are no new user sequences to add." << endl;
	    ufafile.close();
	}
	cout << "total size of user updated:" << user_seqs->size() << endl;
	for(int j=0;j<user_seqs->size();j++){
	    cout << "adding: " << user_seqs->at(j).get_id() << endl;
	}
	gene_db.add_user_seqs_to_db(user_seqs);
	if(user_seqs->size() + keep_seqs->size() == 0)
	  exit(0);
	//if update is more than 20% then redo run
	if(user_seqs->size() + keep_seqs->size() >= (stored_seqs.size() * 0.2)){
	    updateDB = false;
	    cout << "update more than 20% of original size, redoing the whole thing" << endl;
	    return 2;
	}
	//add the right gi numbers before add the rest of the seqs are added to keep_seqs
	write_user_numbers();
	ufafile.close();
	//check files for existing taxonomic break down as it will generally be 
	//these seperations or more fine
	vector<string> align_names;
	gene_db.get_alignment_names(align_names);
	if(usertree == false){
	    for (unsigned int i = 0;i < align_names.size();i++){
		string sql = "SELECT ncbi_id,left_value,right_value FROM taxonomy where ncbi_id = "+align_names[i]+";";
		Query query(conn);
		query.get_result(sql);
		string t_id; string l_id; string r_id;
		while(query.fetch_row()){
		    t_id = query.getstr();
		    l_id = query.getstr();
		    r_id = query.getstr();
		}
		if (t_id.size() == 0){
		    cout << "There is an error getting the id for " << align_names[i] << endl;
		    exit(0);
		}
		for(unsigned int j = 0; j < keep_seqs->size(); j++){
		    //start here, need to get the higher taxon
		    sql = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = '"+keep_seqs->at(j).get_ncbi_tax_id()+"';";
		    Query query2(conn);
		    query2.get_result(sql);
		    string lefttid;
		    string righttid;
		    while(query2.fetch_row()){
			lefttid = query2.getstr();
			righttid = query2.getstr();
		    }
		    if (lefttid > l_id && righttid < r_id){
			if ((int) count(sname_ids.begin(),sname_ids.end(),t_id) == 0){
			    sname_ids.push_back(t_id);
			    snames.push_back(align_names[i]);
			    //add the sequences from the file into keep_seqs , this should be easier when moving to sqlite
			    add_seqs_from_db_to_seqs_vector(align_names[i],keep_seqs,stored_seqs);
			    gene_db.remove_alignment_by_name(align_names[i]);
			    cout << "deleting " << align_names[i] << endl;
			}
			break;
		    }
		}
	    }
	    cout <<"seqs to process: " << keep_seqs->size() << endl;
	}else{//usertree == true
	    //can do seq similarity, tree building distance, or ncbi
	    //going to try ncbi first
	    //basically going to ask what are the ncbi taxa in each, if parent id of any of these is the same
	    //going to test seq similarity and take best, if fails, then will just add to singletons

	    //need to check if the tree is the same as the previous
	    //if the nodes are then I can rename the new tree nodes to those
	    vector<string> rem_files;
	    //basic idea is to take the whole set of sequence names and see if the mrca of the names in a file include any seqs outside of that clade
	    //storedseqs are the seqs from before
	    FastaUtil fu;
	    map<Node *,string> rename_nodes;//this should allow for renaming
            //the file_names left over are all the ones that need to be redone
	    for(int j=0;j<align_names.size();j++){
		bool test = true;
		vector<Sequence> tempseqs;
		fu.readFile(gene_name+"/"+align_names[j],tempseqs);
		//TODO: if one seq need to just delete
		vector<string> mrca_names;
		//need to correct the names for user or not
		for(int i=0;i<tempseqs.size();i++){
		    mrca_names.push_back(tempseqs[i].get_id());
		}
		Node * tmrca = userguidetree->getMRCA(mrca_names);
		if(tmrca == NULL){
		    test = false;
		}else{
		    vector<Node *> leaves=tmrca->get_leaves();
		    for(int i=0;i<leaves.size();i++){
			if((int)count(mrca_names.begin(),mrca_names.end(),leaves[i]->getComment())==0 && 
			   (int)count(stored_seqs_ncbis.begin(),stored_seqs_ncbis.end(),leaves[i]->getComment())>0){
			    test = false;
			    break;
			}
		    }
		}
		if (test == true){
		    tmrca->setName(align_names[j]);
		    rename_nodes[tmrca] = align_names[j];
		}else{
		  add_seqs_from_db_to_seqs_vector(align_names[j],keep_seqs,stored_seqs);
		  gene_db.remove_alignment_by_name(align_names[j]);
		  rem_files.push_back(align_names[j]);
		}
	    }
	    //remove 
	    for(int i=0;i<rem_files.size();i++){
		vector<string>::iterator it;
		it = find(align_names.begin(), align_names.end(),rem_files[i]);
		//align_names.erase(it);
	    }
	    if (align_names.size()>0){//explode is just a redo
		for(unsigned int i=0;i<keep_seqs->size();i++){
		    int bestind;
		    int bestscore=0;
		    for(unsigned int j=0;j<align_names.size();j++){
			vector<Sequence> tempseqs;
			fu.readFile(align_names[j],tempseqs);
			int tscore = get_single_to_group_seq_score(keep_seqs->at(i),tempseqs);
			if (tscore > bestscore)
			    bestind = j;
		    }
		    if ((int) count(snames.begin(),snames.end(),align_names[bestind]) == 0){
			snames.push_back(align_names[bestind]); //need to redo this file
			add_seqs_from_db_to_seqs_vector(align_names[bestind],keep_seqs,stored_seqs);
			gene_db.remove_alignment_by_name(align_names[bestind]);
		    }
		}
	    }else{
		//rename the nodes with the remaped , the others can be random
		int ccount =0;
		for(int i=0;i<userguidetree->getInternalNodeCount();i++){
		    if(rename_nodes.count(userguidetree->getInternalNode(i))==1){
			userguidetree->getInternalNode(i)->setName(rename_nodes[userguidetree->getInternalNode(i)]);
		    }else{
			while((int)count(align_names.begin(),align_names.end(),to_string(ccount))>0){
			    ccount += 1;
			}
			userguidetree->getInternalNode(i)->setName(to_string(ccount));
		    }
		}
		//setup for restart
		snames.push_back(userguidetree->getRoot()->getName());
	    }
	}
    }

    //saturation tests
    if (updateDB == true) {
	saturation_tests(sname_ids, snames, keep_seqs);
    }else{
	write_gi_numbers(keep_seqs);
	gifile.close();
	gene_db.add_seqs_to_db(keep_seqs);
	if(userfasta == true){
	  write_user_numbers();
	  ufafile.close();
	  gene_db.add_user_seqs_to_db(user_seqs);
	}
	sname_ids.push_back(sname_id);
	snames.push_back(clade_name);
	saturation_tests(sname_ids, snames, keep_seqs);
    }

    logfile.close();
    delete known_seqs;
    delete keep_seqs;
    return 0;
}

string SQLiteConstructor::get_cladename(){
    return clade_name;
}

vector <string> SQLiteConstructor::get_search(){
    return search;
}

string SQLiteConstructor::get_genename(){
    return gene_name;
}

double SQLiteConstructor::get_madcutoff(){
    return mad_cutoff;
}

double SQLiteConstructor::get_coverage(){
    return coverage;
}

double SQLiteConstructor::get_identity(){
    return identity;
}

int SQLiteConstructor::get_numthreads(){
    return numthreads;
}

//private

//
/*
 * should retrieve all the matches for a sequence based on the description
 * should return the vector of two strings of the ids,taxon_ids in the sequence table
 */
void SQLiteConstructor::first_seq_search_for_gene_left_right(vector<vector<string> > & retvals){
    Database conn(db);
    string sql;
    if (search.size() == 1){
	sql = "SELECT id,ncbi_id FROM sequence WHERE description LIKE '%"+search[0]+"%'";
    }else{
	sql = "SELECT id,ncbi_id FROM sequence WHERE ";
	for(int i=0;i<search.size()-1;i++){
	    sql = sql + "description LIKE '%"+search[i]+"%' OR ";
	}
	sql = sql + "description LIKE '%"+search[search.size()-1]+"%';";
    }
    //string sql = "SELECT * FROM bioentry WHERE description LIKE '%"+search+"%'";
//	sql = "SELECT * FROM bioentry WHERE description LIKE '%"+search+"%' OR description LIKE '%trnK%'"
    Query query(conn);
    query.get_result(sql);
    int co = 1;
    while(query.fetch_row()){
	vector<string> vals;
	vals.push_back(to_string(query.getval()));
	vals.push_back(to_string(query.getval()));
	retvals.push_back(vals);
	co += 1;
    }
    cout << "search number " << co << endl;
    query.free_result();
}

vector<Sequence> SQLiteConstructor::first_get_seqs_for_name_use_left_right
									(int name_id, vector<vector<string> > & results){
    Database conn(db);
    string deeptaxid;
    deeptaxid = to_string(name_id);
    string sql = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = "+deeptaxid;
    Query query(conn);
    query.get_result(sql);
    int left_value_name;
    int right_value_name;
    while(query.fetch_row()){
	left_value_name = query.getval();
	right_value_name = query.getval();
	main_left = left_value_name;
	main_right = right_value_name;
    }
    query.free_result();
    vector <Sequence> seqs;
    int left_value, right_value,ncbi_id;
    string bioentid;
    for (int i =0 ; i < results.size(); i++){
	string ncbi = results[i][1];
	string taxid = "";
	string edname;
	sql = "SELECT id,left_value,right_value,edited_name FROM taxonomy WHERE ncbi_id = "+ncbi+" and name_class = 'scientific name';";
	Query query2(conn);
	query2.get_result(sql);
	while(query2.fetch_row()){
	    taxid = query2.getval();
	    left_value = query2.getval();
	    right_value = query2.getval();
	    edname = query2.getstr();
	}
	query2.free_result();
	if (left_value > left_value_name && right_value < right_value_name){
	    bioentid = results[i][0];
	    sql = "SELECT accession_id,identifier,description,seq FROM sequence WHERE id = "+bioentid;
	    Query query3(conn);
	    query3.get_result(sql);
	    string descr,acc,gi,sequ;
	    while(query3.fetch_row()){
		acc = query3.getstr();
		gi = query3.getstr();
		descr = query3.getstr();
		sequ = query3.getstr();
	    }
	    query3.free_result();
	    Sequence tseq(ncbi, sequ);
	    tseq.set_ncbi_tax_id(ncbi);
	    tseq.set_ncbi_gi_id(gi);
	    tseq.set_name(edname);
	    seqs.push_back(tseq);
	}
    }
    return seqs;
}

/*
 * this take the literal name from the file so this should be prefiltered
 * to be something that ncbi will match with the names not the edited names
 */
vector<Sequence> SQLiteConstructor::use_only_names_from_file(vector<Sequence> & seqs){
    Database conn(db);
    vector<string> * taxa =new vector<string>();
    vector<string> * taxa_ids =new vector<string>();
    //read file
    ifstream ifs(listfilename.c_str());
    string line;
    while(getline(ifs,line)){
	TrimSpaces(line);
	taxa->push_back(line);
	string sql = "SELECT ncbi_id FROM taxonomy WHERE name = '"+line+"'";
	Query query1(conn);
	query1.get_result(sql);
	while(query1.fetch_row()){
	    int id;
	    id = query1.getval();
	    taxa_ids->push_back(to_string(id));
	}
	query1.free_result();
    }
    cout << taxa_ids->size() << endl;
    /*
     * the adding of wild taxa (this will add ALL the children of the names in the list file)
     */
    if (containswild == true){
	cout << "this file contains higher taxa" << endl;
	vector<string> new_ids;
	for(int i=0;i<taxa_ids->size();i++){
	    //first get the left and right value for the taxa
	    string sql = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = "+taxa_ids->at(i)+";";
	    Query query2(conn);
	    query2.get_result(sql);
	    string lefts;string rights;
	    while(query2.fetch_row()){
		int left = query2.getval();
		int right = query2.getval();
		lefts = to_string(left);
		rights = to_string(right);
	    }
	    sql = "SELECT ncbi_id FROM taxonomy WHERE left_value > "+lefts+" AND right_value < "+rights+" AND name_class = 'scientific name';";
	    Query query(conn);
	    query.get_result(sql);
	    long count = query.num_rows();
	    //exit(0);
	    if (count == 0){
		continue;
	    }else{
		while(query.fetch_row()){
		    new_ids.push_back(to_string(query.getval()));
		}
	    }
	    query.free_result();
	}
	for(int i=0;i<new_ids.size();i++){
	    taxa_ids->push_back(new_ids[i]);
	}
    }
    cout << taxa_ids->size() << endl;
    /*
     * end of the wild taxa
     */

    cout << taxa_ids->size() << " names in the file" << endl;
    ifs.close();
    //end read file
    vector<Sequence> seqs_fn;
    for(int i=0;i<seqs.size();i++){
	string taxid = seqs[i].get_ncbi_tax_id();
	int scount = count(taxa_ids->begin(),taxa_ids->end(),taxid);
	if(scount > 0){
	    seqs_fn.push_back(seqs[i]);
	}
    }
    /*
     * added for higher taxa
     */
    if(containshigher == true && containswild == false){
	cout << "this file contains higher taxa" << endl;
	for(int i=0;i < taxa_ids->size();i++){
	    string sql = "SELECT ncbi_id FROM taxonomy WHERE parent_ncbi_id = "+taxa_ids->at(i)+" and name_class = 'scientific name';";
	    Query query(conn);
	    query.get_result(sql);
	    long count = query.num_rows();
	    cout << count << endl;
	    if(count == 0){
		continue;
	    }else{
		try{
		    Sequence tse = add_higher_taxa(taxa_ids->at(i),seqs);
		    seqs_fn.push_back(tse);
		}catch(int a){

		}
	    }
	    query.free_result();
	}
    }
    /*
     * end added higher taxa
     */
    delete taxa;
    delete taxa_ids;
    return seqs_fn;
}

/*
 * called from use_only_names_from_file,  includes the procedure
 * to deal with higher taxa -- this does a mini phlawd construct
 * on the higher taxa to pick the best one
 * steps are
 * 1) send higher taxa name from use_only_names_from_file (which sends
 * based on a name in the list having children)
 * 2) get all seqs of the higher taxa from seqs (sent along)
 * 3) ortho check again known files
 * 4) take best seq (best overall)
 * 5) store in highertaxa container
 * 6) add these to those that pass the ortho check later (need to store
 * the seq in a file so that the higher taxa id is associated with the
 * smaller one (higher is in the phlawd file, smaller is store) saved
 * with the seq)
 */

Sequence SQLiteConstructor::add_higher_taxa(string taxon_id,vector<Sequence>& seqs){
    vector<string> children_ids = get_final_children(taxon_id);
    //get all the seqs in the set that are within the higher taxon
    vector<Sequence> seqs_fn2;
    for(int i=0;i<seqs.size();i++){
	string taxid = seqs[i].get_ncbi_tax_id();
	int scount = count(children_ids.begin(),children_ids.end(),taxid);
	if(scount > 0){
	    seqs_fn2.push_back(seqs[i]);
	}
    }
    if (seqs_fn2.size() == 0){
	throw 0;
    }else{
	/*
	 * now get the best sequence in the set
	 */
	vector<Sequence> * keep_seqs2 = new vector<Sequence>();
	get_same_seqs_openmp_SWPS3(seqs_fn2,keep_seqs2);
	//take keep_seqs and the known_seqs and get the distances and get the best
	vector<int> scores;
	SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
	for(int i=0;i<known_seqs->size();i++){
	    scores.push_back(get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(i),&known_seqs->at(i)));
	}

	int bestid = 0;
	double bestiden = 0;
	if(keep_seqs2->size() > 0){
	    for (int i=0;i<keep_seqs2->size();i++){
		Sequence tseq = keep_seqs2->at(i);
		double maxiden = 0;
		bool rc = false;
		for (int j=0;j<known_seqs->size();j++){
		    bool trc = false;
		    int ret = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), & tseq);
		    double tsc = double(ret)/double(scores[j]);
		    Sequence tseqrc;
		    tseqrc.set_id(tseq.get_id());
		    tseqrc.set_sequence(tseq.reverse_complement());
		    int retrc = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), &tseqrc);
		    if(retrc > ret){
			trc = true;
			tsc = double(retrc)/double(scores[j]);
		    }
		    if (tsc > maxiden){
			maxiden = tsc;
		    }
		    //cout << tsc << endl;
		}
		if (maxiden >= bestiden){
		    bestid = i;
		    bestiden = maxiden;
		}
	    }
	    Sequence bestseq = keep_seqs2->at(bestid);
	    bestseq.set_id(taxon_id); //name will be the higher taxon name
	    cout << "higher taxa change" << endl;
	    cout << taxon_id << "=" << bestseq.get_ncbi_tax_id() << endl;
	    logfile << "higher taxa change\n";
	    logfile << "ncbi: " << bestseq.get_ncbi_tax_id() << "\n";
	    logfile << "name: " << bestseq.get_name() << "\n";
	    delete keep_seqs2;

	    /*
	     * return the best sequence
	     */
	    return bestseq;
	}else{
	    throw 0;
	}
    }
}

/*
 * excluding taxa from sequences
 */
vector<Sequence> SQLiteConstructor::exclude_names_from_file(vector<Sequence>& seqs){
    Database conn(db);
    vector<string> * taxa_ids =new vector<string>();
    //read file
    ifstream ifs(excludefilename.c_str());
    string line;
    while(getline(ifs,line)){
	TrimSpaces(line);
	string sql;
	if (line[0]=='*') { //this indicates a wildcard and will ignore any taxa with this in the name
	    string trimline = line.substr(1,line.size());
	    sql = "SELECT ncbi_id FROM taxonomy WHERE left_value > "+int_to_string(main_left)+" AND right_value < "+int_to_string(main_right)+" AND name like '%"+trimline+"%'";
	}else{

	    sql = "SELECT ncbi_id FROM taxonomy WHERE name = '"+line+"'";
	}
	Query query1(conn);
	query1.get_result(sql);
	while(query1.fetch_row()){
	    int id;
	    id = query1.getval();
	    taxa_ids->push_back(to_string(id));
	}
	query1.free_result();
    }
    cout << taxa_ids->size() << " names in the file" << endl;
    ifs.close();
    //end read file
    vector<Sequence> seqs_fn;
    for(int i=0;i<seqs.size();i++){
	string taxid = seqs[i].get_ncbi_tax_id();
	bool found = false;
	for(int j=0;j<taxa_ids->size();j++){
	    if (taxa_ids->at(j) == taxid){
		found = true;
		break;
	    }
	}
	if(found == false){
	    seqs_fn.push_back(seqs[i]);
	}
    }
    delete taxa_ids;
    return seqs_fn;
}

/*
 * excluding gis from sequences
 */
vector<Sequence> SQLiteConstructor::exclude_gis_from_file(vector<Sequence> &seqs){
    vector<string> * gi_ids =new vector<string>();
    //read file
    ifstream ifs(exclude_gi_filename.c_str());
    string line;
    while(getline(ifs,line)){
	TrimSpaces(line);
	gi_ids->push_back(line);
    }
    cout << gi_ids->size() << " gis in the file" << endl;
    ifs.close();
    //end read file
    vector<Sequence> seqs_fn;
    for(int i=0;i<seqs.size();i++){
	string giid = seqs[i].get_ncbi_gi_id();
	bool found = false;
	for(int j=0;j<gi_ids->size();j++){
	    if (gi_ids->at(j) == giid){
		found = true;
		break;
	    }
	}
	if(found == false){
	    seqs_fn.push_back(seqs[i]);
	}
    }
    delete gi_ids;
    return seqs_fn;
}

/*
 * include only gis from file
 */
vector <Sequence> SQLiteConstructor::include_gis_from_file(vector<Sequence> & seqs){
    vector<string> * gi_ids = new vector<string> ();
    //read file
    ifstream ifs(include_gi_filename.c_str());
    string line;
    while(getline(ifs,line)){
	TrimSpaces(line);
	gi_ids->push_back(line);
    }
    cout << gi_ids->size() << " gis in the file" << endl;
    ifs.close();
    //end read file
    vector<Sequence> seqs_fn;
    for(int i=0;i<seqs.size();i++){
	//string giid = seqs[i].get_accession();
	string giid = seqs[i].get_ncbi_gi_id();
	bool found = false;
	for(int j=0;j<gi_ids->size();j++){
	    if(gi_ids->at(j)==giid){
		found = true;
		break;
	    }
	}
	if(found == true){
	    seqs_fn.push_back(seqs[i]);
	}
    }
    delete gi_ids;
    return seqs_fn;
}



/*
 * OPENMP version
 */
void SQLiteConstructor::get_same_seqs_openmp_SWPS3_justquery(vector<Sequence> & seqs,  vector<Sequence> * keep_seqs){
    vector<int> known_scores;
    SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
    //SBMatrix mat = swps3_get_premade_SBMatrix( "EDNAFULL" );
    for(int i=0;i<known_seqs->size();i++){
	known_scores.push_back(get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(i),&known_seqs->at(i)));
    }
    map<Sequence*,bool> keep_seqs_rc_map;
    vector<double> justqueryvec;
    vector<double> justqueryvec2;
#pragma omp parallel for shared(keep_seqs_rc_map,justqueryvec)
    for (int i=0;i<seqs.size();i++){
	double maxide = 0;
	double maxcov = 0;
	bool rc = false;
	for (int j=0;j<known_seqs->size();j++){
	    bool trc = false;
	    int ret = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), &seqs[i]);
	    double tsc = double(ret)/double(known_scores[j]);
	    seqs[i].perm_reverse_complement();//make reverse complement
	    int retrc = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), &seqs[i]);
	    seqs[i].perm_reverse_complement();//make it back to original
	    int setsc = get_swps3_score_and_rc_cstyle(mat,&seqs[i],&seqs[i]);
	    double fsetsc = double(ret)/double(setsc);
//	    cout <<i << " " << j << " " << setsc << " " << ret << " " << retrc << " " << known_scores[j] << " " <<  tsc << " " << double(ret)/double(setsc) << endl;
//	    cout << seqs[i].get_sequence() << endl;
//	    cout << known_seqs->at(j).get_sequence() << endl;
//	    cout << seqs[i].get_sequence().size() << endl;
//	    cout << known_seqs->at(j).get_sequence().size() << endl;
//	    exit(0);
	    if(retrc > ret){
		trc = true;
		tsc = double(retrc)/double(known_scores[j]);
	    }
	    if (tsc > maxide && std::numeric_limits<double>::infinity() != tsc){
		maxide = tsc;
		maxcov = fsetsc;
		rc = trc;
	    }
	    if(maxide >= min(identity+(identity*0.5),.99) && justseqquery == false)
	       break;
	}
	justqueryvec.push_back(maxide);
	justqueryvec2.push_back(maxcov);
	if (maxide >= identity && maxcov >= coverage){
	    keep_seqs_rc_map[&seqs[i]] = rc;
	    //with this we don't have to keep track of rc anymore unless we want to
	    if (rc == true)
		seqs[i].perm_reverse_complement();//the sequence is suppose to be reverse complement
	}
    }
    map<Sequence*,bool>::iterator it;
    for(it = keep_seqs_rc_map.begin(); it != keep_seqs_rc_map.end(); it++){
	keep_seqs->push_back(*(*it).first);
    }
    ofstream outfile;
    outfile.open ((gene_name+".seqquery").c_str(),ios::out);
    for (int i=0;i<justqueryvec.size();i++){
	outfile << justqueryvec[i] << "\t" << justqueryvec2[i] << endl;
    }
    outfile.close();
}

void SQLiteConstructor::get_same_seqs_openmp_SWPS3(vector<Sequence> & seqs,  vector<Sequence> * keep_seqs){
    vector<int> known_scores;
    SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
    //SBMatrix mat = swps3_get_premade_SBMatrix( "EDNAFULL" );
    for(int i=0;i<known_seqs->size();i++){
	known_scores.push_back(get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(i),&known_seqs->at(i)));
    }
    map<Sequence*,bool> keep_seqs_rc_map;
#pragma omp parallel for shared(keep_seqs_rc_map)
    for (int i=0;i<seqs.size();i++){
	double maxide = 0;
	double maxcov = 0;
	bool rc = false;
	for (int j=0;j<known_seqs->size();j++){
	    bool trc = false;
	    int ret = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), &seqs[i]);
	    double tsc = double(ret)/double(known_scores[j]);
	    seqs[i].perm_reverse_complement();//make reverse complement
	    int retrc = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), &seqs[i]);
	    seqs[i].perm_reverse_complement();//make it back to original
	    int setsc = get_swps3_score_and_rc_cstyle(mat,&seqs[i],&seqs[i]);
	    double fsetsc = double(ret)/double(setsc);
//	    cout <<i << " " << j << " " << setsc << " " << ret << " " << retrc << " " << known_scores[j] << " " <<  tsc << " " << double(ret)/double(setsc) << endl;
//	    cout << seqs[i].get_sequence() << endl;
//	    cout << known_seqs->at(j).get_sequence() << endl;
//	    cout << seqs[i].get_sequence().size() << endl;
//	    cout << known_seqs->at(j).get_sequence().size() << endl;
//	    exit(0);
	    if(retrc > ret){
		trc = true;
		tsc = double(retrc)/double(known_scores[j]);
	    }
	    if (tsc > maxide && std::numeric_limits<double>::infinity() != tsc){
		maxide = tsc;
		maxcov = fsetsc;
		rc = trc;
	    }
	    if(maxide >= min(identity+(identity*0.5),.99) && justseqquery == false)
	       break;
	}
	if (maxide >= identity && maxcov >= coverage){
	    keep_seqs_rc_map[&seqs[i]] = rc;
	    //with this we don't have to keep track of rc anymore unless we want to
	    if (rc == true)
		seqs[i].perm_reverse_complement();//the sequence is suppose to be reverse complement
	}
    }
    map<Sequence*,bool>::iterator it;
    for(it = keep_seqs_rc_map.begin(); it != keep_seqs_rc_map.end(); it++){
	keep_seqs->push_back(*(*it).first);
    }
}

void SQLiteConstructor::remove_duplicates_SWPS3(vector<Sequence> * keep_seqs){
    vector<string> ids;
    vector<string> unique_ids;
    int mycount;

    //uses NCBI taxon ids for dups
    
    for(unsigned int i =0; i<keep_seqs->size(); i++){
	ids.push_back(keep_seqs->at(i).get_ncbi_tax_id());
	mycount = 0;
	if(unique_ids.size() > 0){
	    mycount = (int) count (unique_ids.begin(),unique_ids.end(), keep_seqs->at(i).get_ncbi_tax_id());
	}
	if(mycount == 0){
	    unique_ids.push_back(keep_seqs->at(i).get_ncbi_tax_id());
	}
    }
    

    /*
     * get the best score for each known seq
     */
    vector<int> scores;
    SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
    for(int i=0;i<known_seqs->size();i++){
	scores.push_back(get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(i),&known_seqs->at(i)));
    }

    vector<int> remove;
    for(unsigned int i=0;i<unique_ids.size();i++){
	int mycount = 0;
	mycount = (int) count (ids.begin(),ids.end(), unique_ids[i]);
	if(mycount > 1){
	    vector<int> tremove;
	    for (int j=0;j<ids.size();j++){
		if(ids[j] == unique_ids[i]){
		    remove.push_back(j);
		    tremove.push_back(j);
		}
	    }
	    int bestid = 0;
	    double bestiden = 0;
	    for (int j=0;j<tremove.size();j++){
		Sequence tseq = keep_seqs->at(tremove[j]);
		double maxiden = 0;
		for (int j=0;j<known_seqs->size();j++){
		    int ret = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), & tseq);
		    double tsc = double(ret)/double(scores[j]);
		    if (tsc > maxiden){
			maxiden = tsc;
		    }
		}
		if (maxiden >= bestiden){
		    bestid = tremove[j];
		    bestiden = maxiden;
		}

	    }
	    vector<int>::iterator it;
	    it = find(remove.begin(), remove.end(),bestid);
	    remove.erase(it);
	}
    }
    sort(remove.begin(),remove.end());
    reverse(remove.begin(),remove.end());
    for(unsigned int i=0;i<remove.size();i++){
	keep_seqs->erase(keep_seqs->begin()+remove[i]);
    }
}

void SQLiteConstructor::reduce_genomes(vector<Sequence> * keep_seqs){
    /*
     * get the best score for each known seq
     */
    vector<int> scores;
    SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
    for(int j=0;j<known_seqs->size();j++){
	scores.push_back(get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j),&known_seqs->at(j)));
    }
    for(unsigned int i =0; i<keep_seqs->size(); i++){
	if(keep_seqs->at(i).get_sequence().size() > 10000){
	    cout << "shrinking a genome: "<< keep_seqs->at(i).get_id() << endl;
	    Sequence tseq = keep_seqs->at(i);
	    double maxiden = 0;
	    int maxknown = 0;
	    for (int j=0;j<known_seqs->size();j++){
		int ret = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), & tseq);
		double tsc = double(ret)/double(scores[j]);
		if (tsc > maxiden){
		    maxiden = tsc;
		    maxknown = j;
		}
	    }
	    const string tempfile = genefoldername+"genome_shrink";
	    vector<Sequence> sc1; 
	    FastaUtil seqwriter;
	    sc1.push_back(keep_seqs->at(i));
	    //for (int j=0;j<known_seqs->size();j++){
	    sc1.push_back(known_seqs->at(maxknown));
	    //}
	    seqwriter.writeFileFromVector(tempfile,sc1);
	    string cmd = "mafft ";
	    cmd += genefoldername+"genome_shrink > ";
	    cmd += genefoldername+"genome_shrink_aln";
	    cout << cmd << endl;
	    cout << "aligning" << endl;
	    system(cmd.c_str());
	    /*string cmd2 = "phyutility -clean 0.5 -in ";
	    cmd2 += genefoldername+"genome_shrink_aln -out ";
	    cmd2 += genefoldername+"genome_shrink_out";
	    cout << cmd << endl;
	    cout << "cleaning" << endl;
	    system(cmd2.c_str());
	    */
	    //instead of phyutility
	    clean_for_genome();
	    /*
	     * reading in the sequencing and replacing
	     */
	    FastaUtil seqreader;
	    vector<Sequence> sequences;
	    seqreader.readFile((genefoldername+"genome_shrink_out").c_str(), sequences);
	    for (int j=0;j<sequences.size();j++){
		if (sequences.at(j).get_id() ==  keep_seqs->at(i).get_id()){
		    keep_seqs->at(i).set_sequence(sequences.at(j).get_sequence());
		}
	    }
	    cout << "shrunk size: "<< keep_seqs->at(i).get_id() << endl;
	}
    }
}

void SQLiteConstructor::clean_for_genome(){
    double percent = 0.5;//missing more than this, then remove
    FastaUtil fu;
    vector<Sequence> tempalseqs;
    fu.readFile(genefoldername+"genome_shrink_aln",tempalseqs);
    cout << "cleaning seqs" << endl;
    int seqlength = tempalseqs[0].get_sequence().size();
    float fseql = float(tempalseqs.size());
    vector<int> removeem;
    for(int j=0;j<seqlength;j++){
	int gaps = 0;
	for(int i=0;i<tempalseqs.size();i++){
	    if(tempalseqs[i].get_sequence()[j] == '-' || tempalseqs[i].get_sequence()[j] == 'N' || tempalseqs[i].get_sequence()[j] == 'n')
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
		a += tempalseqs[i].get_sequence()[j];
	}
	tempalseqs[i].set_sequence(a);
    }
    remove((genefoldername+"genome_shrink_aln").c_str());
    fu.writeFileFromVector(genefoldername+"genome_shrink_out",tempalseqs);
}

vector<string> SQLiteConstructor::get_final_children(string inname_id){
    //special case where inname_id is root and id is therefore 1
    if (inname_id == "1"){
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
	ids.push_back(inname_id);
	vector<string> keepids;
	Database conn(db);
	//testing if tip
	bool tip = true;
	//end testing
	while(!ids.empty()){
	    string sql = "SELECT ncbi_id FROM taxonomy WHERE parent_ncbi_id = "+ids.back()+" and name_class='scientific name';";
	    ids.pop_back();
	    Query query(conn);
	    query.get_result(sql);
	    //StoreQueryResult R = query.store();
	    if(query.num_rows() > 0){
		tip = false;
	    }
	    while(query.fetch_row()){
		string tid = to_string(query.getval());
		ids.push_back(tid);
		keepids.push_back(tid);
	    }
	    query.free_result();
	}
	if (tip == true){
	    keepids.push_back(inname_id);
	}

	vector<string> allids;
	vector<string>allnames;

	for(int i=0;i<keepids.size();i++){
	    string sql = "SELECT name,name_class FROM taxonomy WHERE ncbi_id = ";
	    sql += keepids[i];
	    Query query(conn);
	    query.get_result(sql);
	    //StoreQueryResult R = query.store();
	    while(query.fetch_row()){
		//string tid = R[j][0].c_str();
		string tn = query.getstr();
		string cln = query.getstr();
		if(cln.find("scientific")!=string::npos && tn.find("environmental")==string::npos && cln.find("environmental")==string::npos){
		    allids.push_back(keepids[i]); //was taxon id, now ncbi id
		    allnames.push_back(tn);
		}
	    }
	    query.free_result();
	}
	return allids;
    }
}

/*
 * getting the final ncbi_ids for the children originating at a node
 */
vector<string> SQLiteConstructor::get_final_children_node(Node * node){
    vector<string> allids;
    vector<Node *> leaves = node->get_leaves();
    for (int i=0;i<leaves.size();i++){
	allids.push_back(leaves[i]->getComment());
    }
    return allids;
}

/*
 * same as get_final_children_node but added the ability to do hierarchical (so the tips in tree may or may not be species)
 */
vector<string> SQLiteConstructor::get_final_children_node_hier(Node * node){
    vector<string> allids;
    vector<Node *> leaves = node->get_leaves();
    for(int i=0;i<leaves.size();i++){
	vector<string> tempids = get_final_children(leaves[i]->getComment());
	for(int j=0;j<tempids.size();j++){allids.push_back(tempids[j]);};
    }
    return allids;
}

void SQLiteConstructor::get_seqs_for_names(string inname_id, vector<Sequence> * seqs, vector<Sequence> * temp_seqs){
    vector<string> final_ids;
    final_ids = get_final_children(inname_id);
    for(unsigned int i=0;i<seqs->size();i++){
	string tid = seqs->at(i).get_ncbi_tax_id();
	int mycount = 0;
	mycount = (int) count (final_ids.begin(),final_ids.end(), tid);
	if(mycount > 0){
	    temp_seqs->push_back(seqs->at(i));
	}
    }
}

void SQLiteConstructor::get_seqs_for_names_user(string inname_id, vector<Sequence> * temp_seqs){
    vector<string> final_ids;
    final_ids = get_final_children(inname_id);
    for(unsigned int i=0;i<user_seqs->size();i++){
	string tid = user_seqs->at(i).get_ncbi_tax_id();//was comment
	int mycount = 0;
	mycount = (int) count (final_ids.begin(), final_ids.end(), tid);
	if (mycount>0){
	    temp_seqs->push_back(user_seqs->at(i));
	}
    }
}

/* for userguidetree
 * this is intended to retrieve all the seqs that are contained below a node
 */
void SQLiteConstructor::get_seqs_for_nodes(Node * node, vector<Sequence> * seqs, vector<Sequence> * temp_seqs){
    vector<string> final_ids;
    //final_ids = get_final_children_node(node);//TODO: see below
    final_ids = get_final_children_node_hier(node);//TODO: decide between these two
    for(unsigned int i=0;i<seqs->size();i++){
	string tid = seqs->at(i).get_ncbi_tax_id();
	int mycount = 0;
	mycount = (int) count (final_ids.begin(),final_ids.end(), tid);
	if(mycount > 0){
	    temp_seqs->push_back(seqs->at(i));
	}
    }
}

void SQLiteConstructor::get_seqs_for_user_nodes(Node * node,  vector<Sequence> * temp_seqs){
    vector<string> final_ids;
    vector<Node *> leaves = node->get_leaves();
    for(unsigned int i=0;i<user_seqs->size();i++){
	if (user_fasta_node_map.count(&user_seqs->at(i)) > 0){
	    int mycount =0;
	    mycount = (int) count (leaves.begin(),leaves.end(),user_fasta_node_map[&user_seqs->at(i)]);
	    if (mycount > 0)
		temp_seqs->push_back(user_seqs->at(i));
	}
    }
}

void SQLiteConstructor::make_mafft_multiple_alignment(vector<Sequence> * inseqs, vector<Sequence> * inuserseqs){
    //make file
    vector<double> retvalues;
    FastaUtil seqwriter1;
    vector<Sequence> sc1;
    for(unsigned int i=0;i<inseqs->size();i++){
	sc1.push_back(inseqs->at(i));
    }
    //TODO: do the user sequence names need to be changed
    for(unsigned int i=0;i<inuserseqs->size();i++){
	sc1.push_back(inuserseqs->at(i));
    }
    const string fn1 = genefoldername+"tempfile";
    seqwriter1.writeFileFromVector(fn1,sc1);

    //make alignment
    string cmd = "mafft ";//--thread " + to_string(omp_get_max_threads());
    cmd += genefoldername+"tempfile > ";
    cmd += genefoldername+"outfile 2> ";
    cmd += genefoldername+"mafft.out";
    cout << "aligning" << endl;
    cout << cmd << endl;
    /*
    FILE *fp = popen(cmd.c_str(), "r" );
    char buff[1000];
    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	string line(buff);
    }
    pclose( fp );
    */
    system(cmd.c_str());
}

double SQLiteConstructor::calculate_MAD_quicktree(){
    string line;
    ofstream outfile;
    outfile.open ((genefoldername+"outfile.stoc").c_str(),ios::out);
    FastaUtil fu;
    vector<Sequence> tempseqs;
    fu.readFile(genefoldername+"outfile",tempseqs);
    for(int i=0;i<tempseqs.size();i++){
	outfile << tempseqs[i].get_id() << "\t" << tempseqs[i].get_sequence() << endl;
    }
    outfile.close();

    string cmd = "quicktree -in a -out m ";
    cmd += genefoldername+"outfile.stoc > ";
    cmd += genefoldername+"dist";
    cout << "calculating distance" << endl;
    system(cmd.c_str());

    vector<double> p_values;
    vector<double> jc_values;

    /*
     * read the matrix
     */
    //string line;
    ifstream pfile ((genefoldername+"dist").c_str());
    vector<string> tokens;
    int nspecies = 0;
    int curspecies = 0;
    bool begin = true;
    if (pfile.is_open()){
	while (! pfile.eof() ){
	    getline (pfile,line);
	    string del("\t ");
	    tokens.clear();
	    Tokenize(line, tokens, del);
	    if(tokens.size() > 1){
		double n1;
		for(int j = curspecies; j < nspecies-1;j++){
		    n1 = atof(tokens.at(j+2).c_str());
		    p_values.push_back(n1);
		    //jc will be (-(3./4.)*math.log(1-(4./3.)*p))
		    jc_values.push_back((-(3./4.)*log(1-(4./3.)*n1)));
		}
		curspecies += 1;
	    }else if (begin == true){
		begin = false;
		TrimSpaces(line);
		nspecies = atoi(line.c_str());
	    }
	}
	pfile.close();
    }
    vector<double>all_abs;
    double med = 0;
    for (unsigned int i=0;i<p_values.size();i++){
	all_abs.push_back(fabs(jc_values[i]-p_values[i]));
    }
    med = median(all_abs);
    cout << "median: " << med << endl;
    vector<double> all_meds;
    for (unsigned int i=0;i<p_values.size();i++){
	all_meds.push_back(fabs(med-all_abs[i]));
    }
    return 1.4826*(median(all_meds));
}



double SQLiteConstructor::calculate_MAD_quicktree_sample(vector<Sequence> * inseqs, vector<Sequence> * inuserseqs){
    srand ( time(NULL) );
    vector<int> rands;
    vector<Sequence> tseqs;
    vector<Sequence> tseqs2;
    if(inseqs->size()>0){
	for(int i=0;i<RANDNUM;i++){
	    int n = rand() % inseqs->size();
	    bool x = false;
	    for(int j=0;j<rands.size();j++){
		if(n == rands[j]){
		    x = true;
		}continue;
	    }
	    if(x == true){
		i--;
	    }else{
		rands.push_back(n);
	    }
	}
	sort(rands.begin(),rands.end());
	for(int i=0;i<RANDNUM;i++){
	    tseqs.push_back(inseqs->at(rands[i]));
	}
    }
    if(inuserseqs->size()>0){
	rands.clear();
	for(int i=0;i<RANDNUM;i++){
	    int n = rand() % inuserseqs->size();
	    bool x = false;
	    for(int j=0;j<rands.size();j++){
		if(n == rands[j]){
		    x = true;
		}continue;
	    }
	    if(x == true){
		i--;
	    }else{
		rands.push_back(n);
	    }
	}
	sort(rands.begin(),rands.end());
	for(int i=0;i<RANDNUM;i++){
	    tseqs2.push_back(inuserseqs->at(rands[i]));
	}
    }
    make_mafft_multiple_alignment(&tseqs,&tseqs2);
    return calculate_MAD_quicktree();
}

/*
 * changed this to accept a set of name_ids and names
 * this should allow for more flexible input when updating
 * species
 *
 * the standard run will simply put one name and one id 
 * in the vectors
 */
void SQLiteConstructor::saturation_tests(vector<string> name_ids, vector<string> names, 
	vector<Sequence> * keep_seqs){

    cout << "starting saturation" << endl;
    vector<Sequence> allseqs; 
    for(int i=0;i<keep_seqs->size();i++){
	allseqs.push_back(keep_seqs->at(i));
    }
    //add user seqs now
    if (userfasta == true){
	for(int i=0;i<user_seqs->size();i++){
	    allseqs.push_back(user_seqs->at(i));
	}
    }
    vector<int> exist_alignments;
    
    vector<string> seq_set_filenames;//here they will be not node names but numbers

    if(ncbi_saturation == true){
	string name;
	string name_id;
	while(!names.empty() && !allseqs.empty()){
	    name_id = name_ids.back();
	    name_ids.pop_back();
	    name = names.back();
	    names.pop_back();
	    vector<Sequence> * temp_seqs = new vector<Sequence>();
	    vector<Sequence> * temp_user_seqs = new vector<Sequence>();
	    temp_seqs->empty();
	    temp_user_seqs->empty();
	    get_seqs_for_names(name_id,keep_seqs,temp_seqs);
	    get_seqs_for_names_user(name_id,temp_user_seqs);
	    if(temp_seqs->size() + temp_user_seqs->size() == 1){
		/*
		 * use to make an orphan but rather just make a singleton file
		 */
		cout << name << ": " << name_id << " " << temp_seqs->size() <<" " << temp_user_seqs->size()<<  " " << allseqs.size() << endl;
		//make file
		for(int i=0;i<temp_seqs->size();i++){
		    temp_seqs->at(i).set_aligned_seq(temp_seqs->at(i).get_sequence());
		    remove_seq_from_seq_vector(&allseqs,temp_seqs->at(i).get_id());
		}
		//user seqs
		for(int i=0;i<temp_user_seqs->size();i++){
		    temp_user_seqs->at(i).set_aligned_seq(temp_user_seqs->at(i).get_sequence());
		    remove_seq_from_seq_vector(&allseqs,temp_user_seqs->at(i).get_id());
		}
		int alignid = gene_db.add_alignment(name_id,temp_seqs, temp_user_seqs);
		if(updateDB==true)
		    gene_db.toggle_alignment_update(alignid);
		exist_alignments.push_back(alignid);
	    }else if (temp_seqs->size() + temp_user_seqs->size()== 0){
		continue;
	    }else{
		/*
		 * multiple sequences
		 */
		cout << name << ": " << name_id << " " << temp_seqs->size() << " "<< temp_user_seqs->size() << endl;
		double mad;
		if(temp_seqs->size() + temp_user_seqs->size() > 2){
		    if (temp_seqs->size() +temp_user_seqs->size() < ALIGNLIMIT){
			make_mafft_multiple_alignment(temp_seqs,temp_user_seqs);
			mad = calculate_MAD_quicktree();
/*		    }else if(temp_seqs->size() +temp_user_seqs->size() < 10000){
			mad = 0;
			for (int i=0;i<10;i++)
			    mad = mad + (calculate_MAD_quicktree_sample(temp_seqs,temp_user_seqs)/10.0);
			mad = mad * 2; //make sure is conservative*/
		    }else{//if it is really big
			mad = mad_cutoff + 1;//make sure it gets broken up
		    }
		}else{
		    make_mafft_multiple_alignment(temp_seqs,temp_user_seqs);
		    mad = 0;
		}
		cout << "mad: "<<mad << endl;
		//if mad scores are good, store result
		if (mad <= mad_cutoff){
		    match_aligned_file(temp_seqs,temp_user_seqs);
		    for(int i=0;i<temp_seqs->size();i++){
			remove_seq_from_seq_vector(&allseqs,temp_seqs->at(i).get_id());
		    }
		    //user seqs
		    for(int i=0;i<temp_user_seqs->size();i++){
			remove_seq_from_seq_vector(&allseqs,temp_user_seqs->at(i).get_id());
		    }
		    int alignid = gene_db.add_alignment(name_id,temp_seqs, temp_user_seqs);
		    if(updateDB==true)
			gene_db.toggle_alignment_update(alignid);
		    exist_alignments.push_back(alignid);
		}
		//if mad scores are bad push the children into names
		else{
		    vector<string>child_ids;
		    Database conn(db);
		    string sql = "SELECT ncbi_id FROM taxonomy WHERE parent_ncbi_id = ";
		    sql += name_id;
		    sql += " and name_class = 'scientific name';";
		    Query query(conn);
		    query.get_result(sql);
		    //StoreQueryResult R = query.store();
		    while(query.fetch_row()){
			string sql2 = "SELECT name,name_class FROM taxonomy WHERE ncbi_id = ";
			string resstr = to_string(query.getval());
			sql2 += resstr;
			sql2 += " and name_class = 'scientific name';";
			Query query2(conn);
			query2.get_result(sql2);
			//StoreQueryResult R2 = query2.store();
			while(query2.fetch_row()){
			    string tn = query2.getstr();
			    string cln = query2.getstr();
			    if(cln.find("scientific")!=string::npos && tn.find("environmental")==string::npos && cln.find("environmental")==string::npos){
				string tid = resstr;
				name_ids.push_back(tid);
				names.push_back(tn);
			    }
			}
			query2.free_result();
		    }
		    query.free_result();
		}
	    }
	    delete (temp_seqs);
	    delete (temp_user_seqs);
	}//END NCBI SATURATION
    }else{//user guide tree
	/*
	 * The idea here is to use the tree structure as the guide for the alignment and the 
	 * breaking up of the groups. So the steps are
	 *
	 * 1) the stack is nodes
	 * 2) pop a node off and get all the seqs that are in that 
	 * 3) do everything, if there is saturation push the children in there
	 * 3a) as part of the do everything bit, the mafft alignments should be recieving the node
	 *     as a guide tree for the alignment
	 */
	//for now ignoring the names that are sent because it is just the root, for update it should be the nodes
	cout << "using user guide tree"<<endl;
	vector<Node *> nodes;
	if (updateDB == true){
	    for(int i=0;i<names.size();i++){
		bool found = false;
		for(int j=0;j<userguidetree->getNodeCount();j++){
		    if (userguidetree->getNode(j)->getName()==names[i]){
			nodes.push_back(userguidetree->getNode(j));
			found = true;
			break;
		    }
		}
		if (found == false){
		    cerr << "problem updating and finding this node : "<<names[i]<<endl;
		    exit(0);
		}
	    }
	}else{
	    nodes.push_back(userguidetree->getRoot());//this is different for update
	}
	while(!nodes.empty()){
	    Node * curnode = nodes.back();
	    nodes.pop_back();
	    vector<Sequence> * temp_seqs = new vector<Sequence>();
	    vector<Sequence> * temp_user_seqs = new vector<Sequence>();
	    temp_user_seqs->empty();
	    temp_seqs->empty();
	    get_seqs_for_nodes(curnode,keep_seqs,temp_seqs);
	    get_seqs_for_user_nodes(curnode,temp_user_seqs);
	    cout <<temp_seqs->size()<< " "<< temp_user_seqs->size() << endl;
	    if(temp_seqs->size()+temp_user_seqs->size() == 1){//just one sequence in the group
		cout << curnode->getName() << " " << temp_seqs->size() << " " << temp_user_seqs->size() << endl;
		//make file
		for(int i=0;i<temp_seqs->size();i++){
		    temp_seqs->at(i).set_aligned_seq(temp_seqs->at(i).get_sequence());
		    remove_seq_from_seq_vector(&allseqs,temp_seqs->at(i).get_id());
		}
		//user seqs
		for(int i=0;i<temp_user_seqs->size();i++){
		    temp_user_seqs->at(i).set_aligned_seq(temp_user_seqs->at(i).get_sequence());
		    remove_seq_from_seq_vector(&allseqs,temp_user_seqs->at(i).get_id());
		}
		int alignid = gene_db.add_alignment(curnode->getName(),temp_seqs, temp_user_seqs);
		if(updateDB==true)
		    gene_db.toggle_alignment_update(alignid);
		exist_alignments.push_back(alignid);
	    }else if (temp_seqs->size() + temp_user_seqs->size() == 0){
		continue;
	    }else{
		/*
		 * multiple sequences
		 */
		cout << curnode->getName() << " " << temp_seqs->size() << endl;
		double mad;
		if(temp_seqs->size()+temp_user_seqs->size() > 2){
		    if (temp_seqs->size() +temp_user_seqs->size()< ALIGNLIMIT){
			//TODO: add input tree for mafft
			make_mafft_multiple_alignment(temp_seqs,temp_user_seqs);
			mad = calculate_MAD_quicktree();
			/*}else if(temp_seqs->size() +temp_user_seqs->size()< 10000){
			//need to make this happen 10 tens and average
			mad = 0;
			for (int i=0;i<10;i++)
			    mad = mad + (calculate_MAD_quicktree_sample(temp_seqs,temp_user_seqs)/10.0);
			mad = mad * 2; //make sure is conservative*/
		    }else{//if it is really big
			mad = mad_cutoff + 1;//make sure it gets broken up
		    }
		}else{
		    make_mafft_multiple_alignment(temp_seqs,temp_user_seqs);
		    mad = 0;
		}
		cout << "mad: "<<mad << endl;
		//if mad scores are good, store result
		if (mad <= mad_cutoff){
		    match_aligned_file(temp_seqs,temp_user_seqs);
		    for(int i=0;i<temp_seqs->size();i++){
			remove_seq_from_seq_vector(&allseqs,temp_seqs->at(i).get_id());
		    }
		    for(int i=0;i<temp_user_seqs->size();i++){
			remove_seq_from_seq_vector(&allseqs,temp_user_seqs->at(i).get_id());
		    }
		    int alignid = gene_db.add_alignment(curnode->getName(),temp_seqs, temp_user_seqs);
		    if(updateDB==true)
			gene_db.toggle_alignment_update(alignid);
		    exist_alignments.push_back(alignid);
		}
		//if mad scores are bad push the children into names
		else{
		    //need to get the children
		    for(int i=0;i<curnode->getChildCount();i++){
			nodes.push_back(curnode->getChild(i));
		    }
		}
	    }
	    delete (temp_seqs);
	    delete (temp_user_seqs);
	}
    }
    /*
     * deal with the singletons
     *
     * singletons should be sequences that either don't have any data in the 
     * tree that is input or in the ncbi database if that is the tree to be 
     * used 
     */
    cout << "leftovers: " << allseqs.size() << endl;

    /*
     * if NCBI taxa are all that is wanted, and they are wanted to be 
     * clean
     */
    cout << "picking where leftovers should go" << endl;
    for (int i=0;i<allseqs.size();i++){
	string name;
	cout <<"adding leftover: "<<allseqs.at(i).get_id()<<endl;
	int bestscore = 0;
	int bestind = 0;
	for(int j=0;j<exist_alignments.size();j++){
	    vector<Sequence> tempseqs;
	    gene_db.get_align_seqs_unaligned(exist_alignments[j],tempseqs);
	    int tscore = get_single_to_group_seq_score(allseqs[i],tempseqs);
	    if (tscore > bestscore){
		bestind = j;
	    }
	}
	//TODO: START here
	//align the new one
	vector<Sequence> tempseqs;
	gene_db.get_align_seqs_unaligned(exist_alignments[bestind],tempseqs);
	Sequence tseq(to_string(allseqs[i].get_sqlite_id()),allseqs[i].get_sequence());
	tempseqs.push_back(tseq);
	vector<Sequence> emptys;
	make_mafft_multiple_alignment(&tempseqs,&emptys);
	match_an_aligned_seq(&allseqs[i]);
	//add then update the aligned seqs
	gene_db.add_seq_to_alignment(exist_alignments[bestind],allseqs[i]);
	//read mafft
	vector<Sequence> alseqs;
	get_aligned_file(&alseqs);
	gene_db.update_align_seqs(exist_alignments[bestind],alseqs);
	if(updateDB==true)
	    gene_db.toggle_alignment_update(exist_alignments[bestind]);
    }
    cout << "finished with sequence processing" << endl;
}

/*
 * comparing a sequences to a group of sequences and returning the best score
 * 
 * this can be helpful when deciding with which group a sequence goes
 *
 * ASSUMPTIONS: the sequences are all pointing the right direction
 */
int SQLiteConstructor::get_single_to_group_seq_score(Sequence & inseq,vector<Sequence> & ginseqs){
    vector<int> scores;
    SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
    for (int i=0;i<ginseqs.size();i++){
	double maxide = 0;
	int ret = swps3_maxscores(mat, &inseq,&ginseqs[i]);
	double tsc = double(ret);
	//cout <<i << " " << j << " " << ret << " " << retrc << " " << known_scores[j] << " " <<  tsc << endl;
	if (std::numeric_limits<double>::infinity() != tsc){
	    scores.push_back(tsc);
	}
    }
    int maxide=0;
    for(int i=0;i<scores.size();i++){
	if(scores[i]>maxide)
	    maxide = scores[i];
    }
}


/*
 * this stores the gi numbers for reference
 */

void SQLiteConstructor::write_gi_numbers(vector<Sequence> * dbs){
    for(int i=0;i<dbs->size();i++){
	//gifile << dbs->at(i).get_tax_id() << "\t"; //don't need this anymore
	gifile << dbs->at(i).get_ncbi_tax_id() << "\t";
	//gifile << dbs->at(i).get_accession() << endl;
	gifile << dbs->at(i).get_ncbi_gi_id() << "\t";
	gifile << dbs->at(i).get_name()<<endl;
    }
}

/*
 * this stores the names for the user seqs, mostly this is for updates
 */
void SQLiteConstructor::write_user_numbers(){
  for(int i=0;i<user_seqs->size();i++){
    ufafile <<user_seqs->at(i).get_name() <<"\t";
    ufafile << user_seqs->at(i).get_ncbi_tax_id()<<"\t";
    ufafile << endl;
  }
}

/*
 * this is primarily used for add the seqs from a file to the dbseq and rc vectors 
 * most for updating alignments
 *
 * all the sequences need to be in the database for this to work
 */
void SQLiteConstructor::add_seqs_from_db_to_seqs_vector(string alignname,vector<Sequence> * keep_seqs, vector<Sequence> & storedseqs){
    FastaUtil fu;
    vector<Sequence> tseqs;
    gene_db.get_align_seq_unaligned_fully_initialized(alignname,tseqs);
    cout << "seqs from " << alignname << ": " << tseqs.size() << endl;
    Database conn(db);
    for(unsigned int i=0;i<tseqs.size();i++){
        size_t found;
	//if it is a user seq
        found = tseqs[i].get_id().find("user_");
        if (found != string::npos){
	    //TODO: not sure if I need this
	    if(usertree == true){
		cout << "matching user fasta seqs to user guide tree" << endl;
		int count = 0;
		for(int i=0;i<userguidetree->getExternalNodeCount();i++){
		    string tname = userguidetree->getExternalNode(i)->getName();
		    cout << tname << endl;
		    if (tname == tseqs[i].get_id() || tname == tseqs[i].get_ncbi_tax_id() || ("user_"+tname == tseqs[i].get_id()) || ("user_"+tname == tseqs[i].get_ncbi_tax_id())){
			user_fasta_node_map[&tseqs[i]] = userguidetree->getExternalNode(i);
			count +=1;
		    }
		}
		cout << "matches: "<<count<< " prop:"<< count/user_seqs->size() << endl;
	    }
	    keep_seqs->push_back(tseqs[i]);
	}else{//it is a dbseq
	    keep_seqs->push_back(tseqs[i]);
	}
    }
}


double SQLiteConstructor::get_usertree_keepseq_overlap(vector<Sequence> * keep_seqs){
    set<string> tree_names;
    for(int i=0;i<userguidetree->getExternalNodeCount();i++){
	tree_names.insert(userguidetree->getExternalNode(i)->getComment());
    }
    double ccount = 0;
    for(int i=0;i<keep_seqs->size();i++){
	//need to change this to be hierachical
	if(tree_names.count(keep_seqs->at(i).get_ncbi_tax_id())==1)
	    ccount+=1;
    }
    return ccount/(double)keep_seqs->size();
}

Tree * SQLiteConstructor::get_user_guide_tree_obj(){
    return userguidetree;
}

bool SQLiteConstructor::get_updatestatus(){
    return updateDB;
}

string SQLiteConstructor::get_genedb(){
    return gene_db_name;
}

void SQLiteConstructor::remove_seq_from_seq_vector(vector<Sequence> * inseqs,string sid){
    int eraseint=0;
    int beforesize = inseqs->size();
    for(int zz=0;zz<inseqs->size();zz++){
	if (sid == (*inseqs)[zz].get_id()){
	    eraseint = zz;
	    break;
	}
    }
    inseqs->erase(inseqs->begin()+eraseint);
    if(inseqs->size() != beforesize-1 ){
	cout << "sequence not erased" << endl;
	exit(0);
    }
}

void SQLiteConstructor::match_aligned_file(vector<Sequence> * temp_seqs, vector<Sequence> * temp_user_seqs){
    FastaUtil fu;
    vector<Sequence> tempalseqs;
    fu.readFile(genefoldername+"outfile",tempalseqs);
    for(int i=0;i<tempalseqs.size();i++){
	bool set = false;
	for(int j=0;j<temp_seqs->size();j++){
	    if(tempalseqs[i].get_id() == temp_seqs->at(j).get_id()){
		temp_seqs->at(j).set_aligned_seq(tempalseqs[i].get_sequence());
		set = true;
		break;
	    }
	}
	if(set == false){
	    for(int j=0;j<temp_user_seqs->size();j++){
		if(tempalseqs[i].get_id() == temp_user_seqs->at(j).get_id()){
		    temp_user_seqs->at(j).set_aligned_seq(tempalseqs[i].get_sequence());
		    set = true;
		    break;
		}
	    }
	}
	if(set == false){
	    cout << "error, aligned seq " << tempalseqs[i].get_id() << " has no match" << endl;
	    exit(0);
	}
    }
}

/*
 * assumes the id from the file is the sqlite id from the seq
 */
void SQLiteConstructor::match_an_aligned_seq(Sequence * temp_seq){
    FastaUtil fu;
    vector<Sequence> tempalseqs;
    fu.readFile(genefoldername+"outfile",tempalseqs);
    bool set = false;
    for(int i=0;i<tempalseqs.size();i++){
	if(tempalseqs[i].get_id() == to_string(temp_seq->get_sqlite_id())){
	    temp_seq->set_aligned_seq(tempalseqs[i].get_sequence());
	    set = true;
	    break;
	}
    }
    if(set == false){
	cout << "error, aligned seq " << temp_seq->get_id() << " has no match" << endl;
	exit(0);
    }
}

void SQLiteConstructor::get_aligned_file(vector<Sequence> * temp_seqs){
    FastaUtil fu;
    vector<Sequence> tempalseqs;
    fu.readFile(genefoldername+"outfile",tempalseqs);
    for (int i=0;i<tempalseqs.size();i++){
	tempalseqs[i].set_sqlite_id(atoi(tempalseqs[i].get_id().c_str()));
	tempalseqs[i].set_aligned_seq(tempalseqs[i].get_sequence());
	temp_seqs->push_back(tempalseqs[i]);
    }
}

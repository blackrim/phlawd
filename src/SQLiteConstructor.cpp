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

using namespace std;

#include "omp.h"

#include "libsqlitewrapped.h"


#include "sequence.h"
#include "fasta_util.h"

#include "SQLiteConstructor.h"
#include "DBSeq.h"
#include "utils.h"

//public

template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

SQLiteConstructor::SQLiteConstructor(string cn, vector <string> searchstr, string genen,
		double mad_cut,double cover, double ident, string dbs, string known_seq_filen, bool its, int numt,bool autom,
		bool inupdatedb,string inupdatefile){
    clade_name=cn;
    search=searchstr;
    gene_name = genen;
    mad_cutoff = mad_cut;
    coverage = cover;
    identity=ident;
    db = dbs;
    useITS = its;
    numthreads = numt;
    automated = autom;
    FastaUtil seqReader;
    //added updating seqs from db
    updateDB = inupdatedb;
    //added updating seqs from file
    if(inupdatefile.length() > 0){
	updateFILE = true;
	updatef = inupdatefile;
    }else{
	updateFILE = false;
    }
    known_seqs = new vector<Sequence>();
    seqReader.readFile(known_seq_filen, *known_seqs);
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

void SQLiteConstructor::run(){
    string logn = gene_name;
    logn.append(".log");
    logfile.open(logn.c_str());
    string gin = gene_name;
    gin.append(".gi");
    //updatedb code
    //need to store the existing sequences if there are any
    //TODO: skip the test below for taxa that are already in the file -- add to the exclude names file maybe
    map<string,string> stored_seqs;
    write_EDNAFILE();//if it doesn't exist
    if(updateDB == true){
	cout << "processing existing seqs" << endl;
	gifile.open(gin.c_str(),fstream::in);
	string line;
	bool first = true;
	while(getline(gifile,line)){
	    if (first == true){
		first = false;
		continue;
	    }
	    vector<string> searchtokens;
	    Tokenize(line,searchtokens, "\t");
	    for(int j=0;j<searchtokens.size();j++){
		TrimSpaces(searchtokens[j]);
	    }
	    stored_seqs[searchtokens[0]] = searchtokens[1];
	}
	cout << "existing seqs: " << stored_seqs.size() << endl;
	gifile.close();
	gifile.open(gin.c_str(),fstream::app | fstream::out);
    }else{
	gifile.open(gin.c_str(),fstream::out);
	//gifile << "tax_id\tncbi_tax_id\tgi_number" << endl;
	gifile << "ncbi_tax_id\tgi_number" << endl;
    }
    //end the gi reading

    // if temp directory doesn't exist
    mkdir("TEMPFILES",S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH);
    mkdir(gene_name.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH);

    //if directory exists, delete all files if not updating
    if(updateDB == false){
	//delete all the files in the directory if they are there
	vector<string> exist_filenames;
	getdir(gene_name.c_str(),exist_filenames);
	cout << "removing existing files" << endl;
	for(unsigned int i=0;i<exist_filenames.size();i++){
	    string tname = gene_name+"/"+exist_filenames[i];
	    cout << "removing: " << tname << endl;
	    remove(tname.c_str());
	}
    }

    vector<vector<string> > start_res;
    first_seq_search_for_gene_left_right(start_res);

    //make connection to database
    Database conn(db);
    cout << "connected to " << db << endl;
    vector<int> R;
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
    string sname_id = to_string(R[0]).c_str();
    cout << "Will be using " << name_id << endl;

    //start with a set of seqs given the first clade name and the regions
    vector<DBSeq> startseqs = first_get_seqs_for_name_use_left_right(name_id, start_res);

    cout << "first: " << startseqs.size() << endl;

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

    //use blast to idenify seqs and rc
    vector<DBSeq> * keep_seqs = new vector<DBSeq>();
    //vector<bool> * keep_rc = new vector<bool>();

    /*
     * not sure where ITS should go, but maybe before blasting
     * then the blasting statistics can be strict
     */

    if(useITS == true){
	cout << "starting ITS mode" <<endl;
	combine_ITS(&startseqs);
	cout << "after ITS " << startseqs.size() << endl;
    }
    
    /*
     * comparing the sequences
     */
    get_same_seqs_openmp_SWPS3(startseqs,keep_seqs);
    cout << "blasted: "<< keep_seqs->size() << endl;

    //remove duplicate names
    remove_duplicates_SWPS3(keep_seqs);
    cout << "dups: "<< keep_seqs->size() << endl;
    /*
     * reduce genome sequences
     */
    reduce_genomes(keep_seqs);
	
    //add the updatedb code here
    //get the list of files and record the higher taxa name and 
    //add the additional sequences to the right hierarchy
    vector<string> sname_ids;
    vector<string> snames;
    if(updateDB == true){
	vector<int> toremove;
	for(int j=0;j<keep_seqs->size();j++){
	    if(stored_seqs.count(keep_seqs->at(j).get_ncbi_taxid()) > 0){
		//cout << "removing: " << keep_seqs->at(j).get_ncbi_taxid() << endl;
		toremove.push_back(j);
	    }
	}
	for(int j=0;j<toremove.size();j++){
	    keep_seqs->erase(keep_seqs->begin()+toremove[toremove.size()-(j+1)]);
	}
	//end the program if there is nothing to update
	if(keep_seqs->size() == 0 ){
	    cout << "There are no new sequences to add." << endl;
	    gifile.close();
	    exit(0);//maybe break out to profile
	}
	cout << "total size of updated:" << keep_seqs->size() << endl;
	for(int j=0;j<keep_seqs->size();j++){
	    cout << "adding: " << keep_seqs->at(j).get_ncbi_taxid() << endl;
	}
	//add the write gi numbers before add the rest of the seqs are added to keep_seqs
	write_gi_numbers(keep_seqs);
	gifile.close();
	//check files for existing taxonomic break down as it will generally be 
	//these seperations or more fine
	vector<string> file_names;
	getdir(gene_name,file_names);
	for (unsigned int i = 0;i < file_names.size();i++){
	    string sql = "SELECT ncbi_id,left_value,right_value FROM taxonomy WHERE name = '"+file_names[i]+"';";
	    Query query(conn);
	    query.get_result(sql);
	    string t_id; string l_id; string r_id;
	    while(query.fetch_row()){
		t_id = query.getstr();
		l_id = query.getstr();
		r_id = query.getstr();
	    }
	    if (t_id.size() == 0){
		cout << "There is an error getting the id for " << file_names[i] << endl;
		exit(0);
	    }
	    for(unsigned int j = 0; j < keep_seqs->size(); j++){
		//start here, need to get the higher taxon
		sql = "SELECT left_value,right_value FROM taxonomy WHERE ncbi_id = '"+keep_seqs->at(j).get_ncbi_taxid()+"';";
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
			snames.push_back(file_names[i]);
			//add the sequences from the file into keep_seqs , this should be easier when moving to sqlite
			add_seqs_from_file_to_dbseqs_vector(file_names[i],keep_seqs,stored_seqs);
			string nfilename = gene_name+"/"+file_names[i];
			remove(nfilename.c_str());
		    }
		    break;
		}
	    }
	}
    }

    //saturation tests
    if (updateDB == true) {
	//using the snames as the list of clade names
	//using the snames_ids as the list of ncbi_taxon_ids
	saturation_tests(sname_ids, snames, keep_seqs);
    }else{
	write_gi_numbers(keep_seqs);
	gifile.close();
	sname_ids.push_back(sname_id);
	snames.push_back(clade_name);
	saturation_tests(sname_ids, snames, keep_seqs);
    }

    logfile.close();
    delete known_seqs;
    delete keep_seqs;
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

vector<DBSeq> SQLiteConstructor::first_get_seqs_for_name_use_left_right
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
    }
    query.free_result();
    vector <DBSeq> seqs;
    int left_value, right_value,ncbi_id;
    string bioentid;
    for (int i =0 ; i < results.size(); i++){
	string ncbi = results[i][1];
	string taxid = "";
	sql = "SELECT id,left_value,right_value FROM taxonomy WHERE ncbi_id = "+ncbi+" and name_class = 'scientific name';";
	Query query2(conn);
	query2.get_result(sql);
	while(query2.fetch_row()){
	    taxid = query2.getval();
	    left_value = query2.getval();
	    right_value = query2.getval();
	}
/*
  StoreQueryResult temp1R = query2.store();
  try{
  left_value = temp1R[0][6];
  right_value = temp1R[0][7];
  }catch(mysqlpp::BadConversion E){
  cout << "some sort of taxa [left/right] error (SQLiteConstructor line:267 -- taxon.taxon_id " << taxid <<")" << endl;
  }
*/
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
	    DBSeq tseq(ncbi, sequ, acc, gi,ncbi, taxid, descr);
	    seqs.push_back(tseq);
	}
    }
    return seqs;
}

vector<DBSeq> SQLiteConstructor::use_only_names_from_file(vector<DBSeq> seqs){
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
    vector<DBSeq> seqs_fn;
    for(int i=0;i<seqs.size();i++){
	//string taxid = seqs[i].get_tax_id();
	string taxid = seqs[i].get_ncbi_taxid();
	int scount = count(taxa_ids->begin(),taxa_ids->end(),taxid);
	if(scount > 0){
	    seqs_fn.push_back(seqs[i]);
	}
    }
    //print len(seqs),len(seqs_fn)
    //for seq in seqs_fn:
    //	seqs.remove(seq)
    //print len(seqs),len(seqs_fn)
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
		    DBSeq tse = add_higher_taxa(taxa_ids->at(i),seqs);
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

DBSeq SQLiteConstructor::add_higher_taxa(string taxon_id,vector<DBSeq> seqs){
    vector<string> children_ids = get_final_children(taxon_id);
    //get all the seqs in the set that are within the higher taxon
    vector<DBSeq> seqs_fn2;
    for(int i=0;i<seqs.size();i++){
	//string taxid = seqs[i].get_tax_id();
	string taxid = seqs[i].get_ncbi_taxid();
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
	vector<DBSeq> * keep_seqs2 = new vector<DBSeq>();
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
		DBSeq tseq = keep_seqs2->at(i);
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
	    DBSeq bestseq = keep_seqs2->at(bestid);
	    bestseq.set_id(taxon_id); //name will be the higher taxon name
	    cout << "higher taxa change" << endl;
	    //cout << taxon_id << "=" << bestseq.get_tax_id() << endl;
	    cout << taxon_id << "=" << bestseq.get_ncbi_taxid() << endl;
	    logfile << "higher taxa change\n";
	    //logfile << taxon_id << "=" << bestseq.get_tax_id() << "\n";
	    logfile << "ncbi: " << bestseq.get_ncbi_taxid() << "\n";
	    logfile << "descr: " << bestseq.get_descr() << "\n";
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
vector<DBSeq> SQLiteConstructor::exclude_names_from_file(vector<DBSeq> seqs){
    Database conn(db);
    vector<string> * taxa =new vector<string>();
    vector<string> * taxa_ids =new vector<string>();
    //read file
    ifstream ifs(excludefilename.c_str());
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
    cout << taxa_ids->size() << " names in the file" << endl;
    ifs.close();
    //end read file
    vector<DBSeq> seqs_fn;
    for(int i=0;i<seqs.size();i++){
	//string taxid = seqs[i].get_tax_id();
	string taxid = seqs[i].get_ncbi_taxid();
	int scount = count(taxa_ids->begin(),taxa_ids->end(),taxid);
	if(scount == 0){
	    seqs_fn.push_back(seqs[i]);
	}
    }
    delete taxa;
    delete taxa_ids;
    return seqs_fn;
}

/*
 * excluding gis from sequences
 */
vector<DBSeq> SQLiteConstructor::exclude_gis_from_file(vector<DBSeq> seqs){
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
    vector<DBSeq> seqs_fn;
    for(int i=0;i<seqs.size();i++){
	//string giid = seqs[i].get_accession();
	string giid = seqs[i].get_gi();
	int scount = count(gi_ids->begin(),gi_ids->end(),giid);
	if(scount == 0){
	    seqs_fn.push_back(seqs[i]);
	}
    }
    delete gi_ids;
    return seqs_fn;
}

/*
 * include only gis from file
 */
vector <DBSeq> SQLiteConstructor::include_gis_from_file(vector<DBSeq> seqs){
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
    vector<DBSeq> seqs_fn;
    for(int i=0;i<seqs.size();i++){
	//string giid = seqs[i].get_accession();
	string giid = seqs[i].get_gi();
	int scount = count(gi_ids->begin(),gi_ids->end(),giid);
	if(scount == 1){
	    seqs_fn.push_back(seqs[i]);
	}
    }
    delete gi_ids;
    return seqs_fn;
}



/*
 * OPENMP version
 */
void SQLiteConstructor::get_same_seqs_openmp_SWPS3(vector<DBSeq> & seqs,  vector<DBSeq> * keep_seqs){
    vector<int> known_scores;
    SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
    //SBMatrix mat = swps3_get_premade_SBMatrix( "EDNAFULL" );
    for(int i=0;i<known_seqs->size();i++){
	known_scores.push_back(get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(i),&known_seqs->at(i)));
    }
    map<DBSeq*,bool> keep_seqs_rc_map;
#pragma omp parallel for shared(keep_seqs_rc_map)
    for (int i=0;i<seqs.size();i++){
	double maxide = 0;
	bool rc = false;
	for (int j=0;j<known_seqs->size();j++){
	    bool trc = false;
	    int ret = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), &seqs[i]);
	    double tsc = double(ret)/double(known_scores[j]);
	    seqs[i].perm_reverse_complement();//make reverse complement
	    int retrc = get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(j), &seqs[i]);
	    seqs[i].perm_reverse_complement();//make it back to original
	    //cout <<i << " " << j << " " << ret << " " << retrc << " " << known_scores[j] << " " <<  tsc << endl;
	    if(retrc > ret){
		trc = true;
		tsc = double(retrc)/double(known_scores[j]);
	    }
	    if (tsc > maxide && std::numeric_limits<double>::infinity() != tsc){
		maxide = tsc;
		rc = trc;
	    }
	}
	if (maxide >= identity){
	    keep_seqs_rc_map[&seqs[i]] = rc;
	    //with this we don't have to keep track of rc anymore unless we want to
	    if (rc == true)
		seqs[i].perm_reverse_complement();//the sequence is suppose to be reverse complement
	}
    }
    map<DBSeq*,bool>::iterator it;
    for(it = keep_seqs_rc_map.begin(); it != keep_seqs_rc_map.end(); it++){
	keep_seqs->push_back(*(*it).first);
    }
}

void SQLiteConstructor::remove_duplicates_SWPS3(vector<DBSeq> * keep_seqs){
    vector<string> ids;
    vector<string> unique_ids;
    int mycount;

    //uses database taxon ids for dups
    for(unsigned int i =0; i<keep_seqs->size(); i++){
	ids.push_back(keep_seqs->at(i).get_id());
	mycount = 0;
	if(unique_ids.size() > 0){
	    mycount = (int) count (unique_ids.begin(),unique_ids.end(), keep_seqs->at(i).get_id());
	}
	if(mycount == 0){
	    unique_ids.push_back(keep_seqs->at(i).get_id());
	}
    }

    //uses NCBI taxon ids for dups
    /*
      for(unsigned int i =0; i<keep_seqs->size(); i++){
      ids.push_back(keep_seqs->at(i).get_ncbi_taxid());
      mycount = 0;
      if(unique_ids.size() > 0){
      mycount = (int) count (unique_ids.begin(),unique_ids.end(), keep_seqs->at(i).get_ncbi_taxid());
      }
      if(mycount == 0){
      unique_ids.push_back(keep_seqs->at(i).get_ncbi_taxid());
      }
      }
    */

    /*
     * get the best score for each known seq
     */
    vector<int> scores;
    SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
    for(int i=0;i<known_seqs->size();i++){
	scores.push_back(get_swps3_score_and_rc_cstyle(mat,&known_seqs->at(i),&known_seqs->at(i)));
    }

    vector<int> remove;
//#pragma omp parallel for shared(remove)
    for(unsigned int i=0;i<unique_ids.size();i++){
	mycount = 0;
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
		DBSeq tseq = keep_seqs->at(tremove[j]);
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

void SQLiteConstructor::reduce_genomes(vector<DBSeq> * keep_seqs){
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
	    DBSeq tseq = keep_seqs->at(i);
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
	    /*
	     * shrink with phyutility
	     */
	    const string tempfile = "TEMPFILES/genome_shrink";
	    vector<Sequence> sc1; 
	    FastaUtil seqwriter;
	    sc1.push_back(keep_seqs->at(i));
	    //for (int j=0;j<known_seqs->size();j++){
	    sc1.push_back(known_seqs->at(maxknown));
	    //}
	    seqwriter.writeFileFromVector(tempfile,sc1);
	    const char * cmd = "mafft --thread 2 --auto TEMPFILES/genome_shrink > TEMPFILES/genome_shrink_aln";
	    cout << "aligning" << endl;
	    FILE *fp = popen(cmd, "r" );
	    char buff[1000];
	    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
		string line(buff);
	    }
	    pclose( fp );
	    const char * cmd2 = "phyutility -clean 0.5 -in TEMPFILES/genome_shrink_aln -out TEMPFILES/genome_shrink_out";
	    cout << "cleaning" << endl;
	    FILE *fp2 = popen(cmd2, "r" );
	    char buff2[1000];
	    while ( fgets( buff2, sizeof buff2, fp2 ) != NULL ) {//doesn't exit out
		string line(buff2);
	    }
	    pclose( fp2 );
	    /*
	     * reading in the sequencing and replacing
	     */
	    FastaUtil seqreader;
	    vector<Sequence> sequences;
	    seqreader.readFile("TEMPFILES/genome_shrink_out", sequences);
	    for (int j=0;j<sequences.size();j++){
		if (sequences.at(j).get_id() ==  keep_seqs->at(i).get_id()){
		    keep_seqs->at(i).set_sequence(sequences.at(j).get_sequence());
		}
	    }
	    cout << "shrunk size: "<< keep_seqs->at(i).get_id() << endl;
	}
    }
}

vector<string> SQLiteConstructor::get_final_children(string inname_id){
    vector<string> ids;
    ids.push_back(inname_id);
    vector<string>	keepids;

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

void SQLiteConstructor::get_seqs_for_names(string inname_id, vector<DBSeq> * seqs, vector<DBSeq> * temp_seqs){
    vector<string> final_ids;
    final_ids = get_final_children(inname_id);
    for(unsigned int i=0;i<seqs->size();i++){
	//string tid = seqs->at(i).get_tax_id();
	string tid = seqs->at(i).get_ncbi_taxid();
	int mycount = 0;
	mycount = (int) count (final_ids.begin(),final_ids.end(), tid);
	if(mycount > 0){
	    temp_seqs->push_back(seqs->at(i));
	}
    }
}

void SQLiteConstructor::make_mafft_multiple_alignment(vector<DBSeq> * inseqs){
    //make file
    vector<double> retvalues;
    FastaUtil seqwriter1;
    vector<Sequence> sc1;
    for(unsigned int i=0;i<inseqs->size();i++){
	sc1.push_back(inseqs->at(i));
    }
    const string fn1 = "TEMPFILES/tempfile";
    seqwriter1.writeFileFromVector(fn1,sc1);

    //make alignment
    string cmd = "mafft --thread " + to_string(omp_get_max_threads());
    cmd += " --auto TEMPFILES/tempfile > TEMPFILES/outfile";
    cout << "aligning" << endl;
    FILE *fp = popen(cmd.c_str(), "r" );
    char buff[1000];
    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	string line(buff);
    }
    pclose( fp );
}
/*
 * should replace paup
 * delete paup when done testing
 */

double SQLiteConstructor::calculate_MAD_quicktree(){
    const char * phcmd = "phyutility -concat -in TEMPFILES/outfile -out TEMPFILES/outfile.nex";
    FILE *phfp = popen(phcmd, "r" );
    pclose( phfp );

    ifstream infile;
    ofstream outfile;
    infile.open ("TEMPFILES/outfile.nex",ios::in);
    outfile.open ("TEMPFILES/outfile.stoc",ios::out);
    bool begin = false;
    bool end = false;
    string line;
    /*
     * convert to stockholm format
     */
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

    vector<double> p_values;
    vector<double> jc_values;

    /*
     * read the matrix
     */
    //string line;
    ifstream pfile ("TEMPFILES/dist");
    vector<string> tokens;
    int nspecies = 0;
    int curspecies = 0;
    begin = true;
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

double SQLiteConstructor::calculate_MAD_quicktree_sample(vector<DBSeq> * inseqs){
	srand ( time(NULL) );
	vector<int> rands;
	for(int i=0;i<1000;i++){
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
	vector<DBSeq> tseqs;
	for(int i=0;i<1000;i++){
		tseqs.push_back(inseqs->at(rands[i]));
	}
	make_mafft_multiple_alignment(&tseqs);
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
	vector<DBSeq> * keep_seqs){

//void SQLiteConstructor::saturation_tests(string name_id, vector<DBSeq> * keep_seqs, vector<bool> * keep_rc){
//	vector<string> name_ids;
//	vector<string> names;
//	name_ids.push_back(name_id);
//	names.push_back(clade_name);

    vector<DBSeq> orphan_seqs;
    vector<Sequence> allseqs; 
    for(int i=0;i<keep_seqs->size();i++){
	allseqs.push_back(keep_seqs->at(i));
    }

    string name;
    string name_id;
    while(!names.empty()){
	name_id = name_ids.back();
	name_ids.pop_back();
	name = names.back();
	names.pop_back();
	vector<DBSeq> * temp_seqs = new vector<DBSeq>();
	temp_seqs->empty();
	get_seqs_for_names(name_id,keep_seqs,temp_seqs);
	if(temp_seqs->size() == 1){
	    /*
	     * use to make an orphan but rather just make a singleton file
	     */
	    cout << name << " " << temp_seqs->size() << endl;
	    //make file
	    FastaUtil seqwriter1;
	    vector<Sequence> sc1;
	    for(int i=0;i<temp_seqs->size();i++){
		//need to implement a better way, but this is it for now
		int eraseint=0;
		for(int zz=0;zz<allseqs.size();zz++){
		    if (temp_seqs->at(i).get_id() == allseqs[zz].get_id()){
			eraseint = zz;
			break;
		    }
		}
		allseqs.erase(allseqs.begin()+eraseint);
		sc1.push_back(temp_seqs->at(i));
	    }
	    string fn1 = gene_name;
	    fn1 += "/" + name;
	    seqwriter1.writeFileFromVector(fn1,sc1);
	}else if (temp_seqs->size() == 0){
	    continue;
	}else{
	    cout << name << " " << temp_seqs->size() << endl;
	    double mad;
	    if(temp_seqs->size() > 2){
		if (temp_seqs->size() < 3000){
		    make_mafft_multiple_alignment(temp_seqs);
		    mad = calculate_MAD_quicktree();
		}else if(temp_seqs->size() < 10000){
		    //need to make this happen 10 tens and average
		    mad = 0;
		    for (int i=0;i<10;i++)
			mad = mad + (calculate_MAD_quicktree_sample(temp_seqs)/10.0);
		    mad = mad * 2; //make sure is conservative
		}else{
		    mad = mad_cutoff + 1;//make sure it gets broken up
		}
	    }else{
		mad = 0;
	    }
	    cout << "mad: "<<mad << endl;
	    //if mad scores are good, store result
	    //write to file and to sqlite database
	    if (mad <= mad_cutoff){
		FastaUtil seqwriter1;
		vector<Sequence> sc1; 
		for(int i=0;i<temp_seqs->size();i++){
		    //need to implement a better way, but this is it for now
		    int eraseint=0;
		    for(int zz=0;zz<allseqs.size();zz++){
			if (temp_seqs->at(i).get_id() == allseqs[zz].get_id()){
			    eraseint = zz;
			    break;
			}
		    }
		    allseqs.erase(allseqs.begin()+eraseint);
		    sc1.push_back(temp_seqs->at(i));
		}
		string fn1 = gene_name;
		fn1 += "/" + name;
		seqwriter1.writeFileFromVector(fn1,sc1);
		//SQLITE database storing

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
    }
    /*
     * deal with the singletons
     */
    cout << "leftovers: " << allseqs.size() << endl;
    for(int i=0;i<allseqs.size();i++){
	Database conn(db);
	vector<Sequence> sc1; 
	sc1.push_back(allseqs.at(i));
	string name;
	string sql = "SELECT name,name_class FROM taxonomy WHERE ncbi_id = ";
	sql += allseqs.at(i).get_id();
	cout <<"-"<<allseqs.at(i).get_id()<<endl;
	Query query(conn);
	query.get_result(sql);
	//StoreQueryResult R = query.store();
	while(query.fetch_row()){
	    string tn = query.getstr();
	    string cln = query.getstr();
	    if(cln.find("scientific")!=string::npos && tn.find("environmental")==string::npos && cln.find("environmental")==string::npos){
		string tid = allseqs.at(i).get_id();
		name = tn;
	    }
	}
	query.free_result();
	cout << name << endl;
	FastaUtil seqwriter1;
	string fn1 = gene_name;
	fn1 += "/" + name;
	seqwriter1.writeFileFromVector(fn1,sc1);
    }
}

/*
 * this stores the gi numbers for reference
 */

void SQLiteConstructor::write_gi_numbers(vector<DBSeq> * dbs){
    for(int i=0;i<dbs->size();i++){
	//gifile << dbs->at(i).get_tax_id() << "\t"; //don't need this anymore
	gifile << dbs->at(i).get_ncbi_taxid() << "\t";
	//gifile << dbs->at(i).get_accession() << endl;
	gifile << dbs->at(i).get_gi() << endl;
    }
}



/*
 * this is primarily used for add the seqs from a file to the dbseq and rc vectors 
 * most for updating alignments
 *
 * all the sequences need to be in the database for this to work
 */
void SQLiteConstructor::add_seqs_from_file_to_dbseqs_vector(string filename,vector<DBSeq> * keep_seqs,map<string,string> & taxgimap){
    FastaUtil fu;
    vector<Sequence> tseqs;
    fu.readFile(gene_name+"/"+filename,tseqs);
    cout << "seqs from " << filename << ": " << tseqs.size() << endl;
    Database conn(db);
    for(unsigned int i=0;i<tseqs.size();i++){
	string ncbi = taxgimap[tseqs[i].get_id()];
	string sql = "SELECT accession_id,description FROM sequence WHERE identifier = "+ncbi+";";
	Query query3(conn);
	query3.get_result(sql);
	string descr,acc;
	while(query3.fetch_row()){
	    acc = query3.getstr();
	    descr = query3.getstr();
	}
	query3.free_result();
	DBSeq tseq(tseqs[i].get_id(), tseqs[i].get_sequence(), acc, ncbi, tseqs[i].get_id(), " ", descr);
	keep_seqs->push_back(tseq);
    }
}



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


using namespace std;

#include "tree.h"
#include "node.h"
#include "tree_utils.h"
#include "tree_reader.h"
#include "sequence.h"
#include "fasta_util.h"


#include "utils.h"

#include "libsqlitewrapped.h"

#include "SQLiteProfiler.h"

template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

SQLiteProfiler::SQLiteProfiler(string gn, string cn, string dbs,bool autom,bool updb){
    gene_name = gn;
    use_orphan = false;
    cladename = cn;
    db = dbs;
    automated = autom;
    updatedb = updb;
    profilefoldername = gene_name+"_PROFILE/";
}

void SQLiteProfiler::prelimalign(){
    // if temp directory doesn't exist
    mkdir(profilefoldername.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IWOTH);
    bool standard = true;
    if(updatedb == true){
	vector<string> samefiles; //files were not updated at all
	bool newfiles=false; // if there are new files from MAD in align stage, kick out and tree as new for now 
	vector<string> originalalnfiles;
	vector<string> curfiles;
	getdir(profilefoldername.c_str(),curfiles);
	getdir(gene_name.c_str(),originalalnfiles);
	//to get updated files, just diff the seq counts for each original file and the profilefiles
	for(unsigned int i=0;i< originalalnfiles.size();i++){
	    int origcounts = count_seqs(gene_name,originalalnfiles[i]);
	    if(count(curfiles.begin(),curfiles.end(),originalalnfiles[i]) == 1){
		int counts = count_seqs(profilefoldername,originalalnfiles[i]);
		if (counts < origcounts){
		    updatedfiles.push_back(originalalnfiles[i]);
		}else if (counts == origcounts){
		    samefiles.push_back(originalalnfiles[i]);
		}
	    }else{//count == 0
		newfiles = true;
		cout << "there is a new file: " << originalalnfiles[i] << endl;
		break;
	    }
	}if (newfiles == true){
	    //delete all the files in the directory if they are there
	    vector<string> exist_filenames;
	    getdir(profilefoldername.c_str(),exist_filenames);
	    cout << "removing existing files" << endl;
	    for(unsigned int i=0;i<exist_filenames.size();i++){
		string tname = profilefoldername+"/"+exist_filenames[i];
		cout << "removing: " << tname << endl;
		remove(tname.c_str());
	    }
	    standard = true;
	    updatedb = false;
	}else{//update the profiles
	    standard = false;
	    //files to skip are record.log, and the original alignment files, except for the updated ones
	    //delete everything else
	    //reading the record.log
	    string tname = profilefoldername+"/record.log";
	    ifstream ifs(tname.c_str());
	    string line;
	    while(getline(ifs,line)){
		vector<string> tokens;
		string del(",");
		tokens.clear();
		Tokenize(line,tokens,del);
		for(int i=0;i<tokens.size();i++){TrimSpaces(tokens[i]);}
		profilerun.push_back(tokens[0]);
		vector<string> tokens2;
		string del2("++");
		tokens2.clear();
		Tokenize(tokens[1],tokens2,del2);
		//look to see if an updated file is in the profile
		for(unsigned int i=0;i<updatedfiles.size();i++){
		    size_t found1 = tokens2[0].find(updatedfiles[i]);
		    size_t found2 = tokens2[1].find(updatedfiles[i]);
		    //cout << found1 << " " << found2 << " " << tokens2[0] << " " << tokens2[1] << " " << updatedfiles[i] << endl;
		    if(found1 != string::npos){// find it in the first part of the profile
			string tname = profilefoldername+"/"+tokens[0];
			cout << "removing profile file: " << tname << endl;
			remove(tname.c_str());
			flaggedprofiles.push_back(tokens[0]);
		    }else if(found2 != string::npos){ //find it in the second part of the profile
			string tname = profilefoldername+"/"+tokens[0];
			cout << "removing profile file: " << tname << endl;
			remove(tname.c_str());
			flaggedprofiles.push_back(tokens[0]);
		    }
		}
		vector<string> tvec;
		tvec.push_back(tokens2[0]);
		tvec.push_back(tokens2[1]);
		profilekey[tokens[0]] =tvec;
	    }
	    //deleting things
	    for(unsigned int i=0;i<updatedfiles.size();i++){
		string tname = profilefoldername+"/"+updatedfiles[i];
		cout << "removing updated file: " << tname << endl;
		remove(tname.c_str());
	    }
	    //delete everything that isn't in the same file or in the profilerun vector
	    for(unsigned int i=0;i<curfiles.size();i++){
		if(count(samefiles.begin(),samefiles.end(),curfiles[i]) == 0){
		    if(count(profilerun.begin(),profilerun.end(),curfiles[i]) == 0){
			if(curfiles[i] != "record.log"){
			    string tname = profilefoldername+"/"+curfiles[i];
			    cout << "removing file: "<< tname << endl;
			    remove(tname.c_str());
			}
		    }
		}
	    }
	}
    }else{//standard run
	//delete all the files in the directory if they are there
	vector<string> exist_filenames;
	getdir(profilefoldername.c_str(),exist_filenames);
	cout << "removing existing files" << endl;
	for(unsigned int i=0;i<exist_filenames.size();i++){
	    string tname = profilefoldername+"/"+exist_filenames[i];
	    cout << "removing: " << tname << endl;
	    remove(tname.c_str());
	}
    }
	
    if(standard == true){//not an update run
	file_names = vector<string>();
	cout << "getting file names" << endl;
	getdir(gene_name,file_names);
	for (unsigned int i = 0;i < file_names.size();i++) {
	    /*
	     * test to see how many sequences in the file
	     * if there is only one, just copy, if there
	     * are many, then align
	     */
	    bool many = false;
	    int intReturn = count_seqs(gene_name, file_names[i]);
	    if(intReturn > 1){
		many = true;
	    }
	    if(many == false){
		string cmd = "cp ";
		cmd += gene_name+"/";
		string tfilen = file_names[i];
		fix_bad_chars(tfilen);
		cmd += tfilen;
		cmd += " ";
		cmd += profilefoldername;
		cmd += tfilen;//fix spaces here
		cout << "copying file with single fasta" << endl;
		cout << cmd << endl;
		FILE *fp = popen(cmd.c_str(), "r" );
		char buff[1000];
		while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
		    string line(buff);
		}
		pclose( fp );
	    }else{
		cout << file_names[i] << endl;
		string cmd = "mafft --thread 2 --auto ";
		cmd += gene_name+"/";
		//cmd += file_names[i];//fix spaces here
		string tfilen = file_names[i];
		fix_bad_chars(tfilen);
		cmd += tfilen;
		cmd += " > ";
		cmd += profilefoldername;
		cmd += tfilen;//fix spaces here
		cout << "prelim aligning" << endl;
		cout << cmd << endl;
		FILE *fp = popen(cmd.c_str(), "r" );
		char buff[1000];
		while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
		    string line(buff);
		}
		pclose( fp );
	    }
	    if(intReturn > 100){
		clean_before_profile(file_names[i]);
	    }
	}
    }else{//updaterun
	//only align those files that are updated
	file_names = vector<string>();
	cout << "aligning updated files" << endl;
	for(unsigned int i=0;i<updatedfiles.size();i++){
	    int intReturn = count_seqs(gene_name,updatedfiles[i]);
	    bool many = false;
	    if(intReturn > 1){
		many = true;
	    }
	    if(many == false){
		string cmd = "cp ";
		cmd += gene_name+"/";
		string tfilen = updatedfiles[i];
		fix_bad_chars(tfilen);
		cmd += tfilen;
		cmd += " ";
		cmd += profilefoldername;
		cmd += tfilen;//fix spaces here
		cout << "copying file with single fasta" << endl;
		cout << cmd << endl;
		FILE *fp = popen(cmd.c_str(), "r" );
		char buff[1000];
		while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
		    string line(buff);
		}
		pclose( fp );
	    }else{
		cout << updatedfiles[i] << endl;
		string cmd = "mafft --thread 2 --auto ";
		cmd += gene_name+"/";
		string tfilen = updatedfiles[i];
		fix_bad_chars(tfilen);
		cmd += tfilen;
		cmd += " > ";
		cmd += profilefoldername;
		cmd += tfilen;//fix spaces here
		cout << "prelim aligning" << endl;
		cout << cmd << endl;
		FILE *fp = popen(cmd.c_str(), "r" );
		char buff[1000];
		while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
		    string line(buff);
		}
		pclose( fp );
	    }
	    if(intReturn > 100){
		clean_before_profile(updatedfiles[i]);
	    }
	}
    }
}

int SQLiteProfiler::count_seqs(string dirc, string file_name){
    string cmd = "grep -c \\> ";
    cmd += dirc+"/";
    string tfilen = file_name;
    fix_bad_chars(tfilen);
    cmd += tfilen;
    FILE *fp = popen(cmd.c_str(), "r" );
    char buff[1000];
    int intReturn;
    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	string line(buff);
	intReturn = atoi(line.c_str());
	break;
    }
    pclose( fp );
    return intReturn;
}

void SQLiteProfiler::run(){
    if(updatedb == false){//standard
	file_names = vector<string>();
	cout << "getting file names" << endl;
	getdir(profilefoldername.c_str(),file_names);
	if(file_names.size() > 1){
	    //get guide tree
	    //if the guide tree is not there, then use the ncbi tree
	    string guiname = gene_name+".guide";
	    bool flag = false;
	    fstream fin;
	    fin.open(guiname.c_str(),ios::in);
	    if( fin.is_open() ){
		flag=true;
	    }
	    fin.close();
	    //end test for guide tree
	    map<string,string> numnames;
	    map<string,string> namesnum;
	    vector< vector<double> > numlist;
	    if(flag == true){
		cout << "user guide tree" << endl;
		Tree * tree = read_user_guide_tree(guiname);
		create_distances(file_names,tree,&numnames,&namesnum,&numlist);
	    }else{//use ncbi tree
		cout << "ncbi guide tree" << endl;
		create_distances(cladename,file_names,&numnames,&namesnum,&numlist);
	    }
	    //start profiling
	    cout<<"profiling"<<endl;
	    profile(numnames,namesnum,numlist);
	}else{
	    //need to create an empty record file
	    string recordname = profilefoldername+"/record.log";
	    ofstream ofs(recordname.c_str());
	    ofs.close();
	    copy_final_file(file_names[0]);
	}
    }else{//updatedb
	update_profile();
    }
    rename_final_alignment("FINAL.aln");//requires FINAL.aln
    remove_outliers();//requires FINAL.aln
    rename_final_alignment("FINAL.aln.cln"); //requires FINAL.aln.cln

}

//for updating alignments
//send  the profile key and a profile strings
//
//return the string if there is no match
//return the profile number (as a string if there is a match)
string SQLiteProfiler::get_profilekey_value(string profile_string){
	bool match = false;
	string match_string;
	for(unsigned int i= 0; i< profilerun.size(); i++){
		vector<string> tstrings = profilekey[profilerun[i]];
		string combined = tstrings[0]+"__"+tstrings[1];
		if (combined == profile_string){
			match = true;
			match_string = profilerun[i];
			break;
		}
	}
	if(match == false){
		return profile_string;
	}
	return match_string;
}

Tree * SQLiteProfiler::read_user_guide_tree(string filen){
    TreeReader nw;
    ifstream infile2(filen.c_str());
    vector<string> lines;
    string line;
    while (getline(infile2, line)){
	lines.push_back(line);
    }
    infile2.close();

    Tree * tree = nw.readTree(lines[0]);
    vector<string> orphans;
    for(int i=0;i<file_names.size();i++){
	try{
	    tree->getExternalNode(file_names[i])->getName();
	}catch(int e){
	    orphans.push_back(file_names[i]);
	    use_orphan = true;
	}
    }
    string orphfilename = profilefoldername+"orphan";
    ofstream myfile(orphfilename.c_str());
    for(int i=0;i<orphans.size();i++){
	ifstream tfile (orphans[i].c_str());
	string line;
	while(getline(tfile,line)){
	    myfile << line; // make sure it adds the \n
	}
    }
    myfile.close();
    return tree;
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
void SQLiteProfiler::create_distances(string clade_name, vector<string> names,map<string,string> * numnames,
		map<string,string>* namesnum, vector< vector<double> > * numlist){
    //the map of clade name and number
    numnames->clear();
    //the opposite
    namesnum->clear();
    //the distances in terms of i and j
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
    for(int i=0;i<names.size();i++){
	Database conn(db);
	sql = "SELECT ncbi_id FROM taxonomy WHERE name = '"+names[i]+"' and name_class = 'scientific name';";
	//cout << sql << endl;
	Query query2(conn);
	query2.get_result(sql);
//		StoreQueryResult R = query2.store();
	string nameid = get_right_one(allids, query2);
	sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = "+nameid+" and name_class = 'scientific name';";
	Query query4(conn);
	query4.get_result(sql);
//		StoreQueryResult R2 = query4.store();
	string parentid = get_right_one(allids,query4);
	vector<string> route;
	while(parentid != cladeid){
	    //cout << "loop" << endl;
	    route.push_back(parentid);
	    cout << "nameid " << nameid << endl;
	    cout << "parentid1 " << parentid << endl;
	    nameid = parentid;
	    sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = "+nameid+" and name_class = 'scientific name';";
	    //cout << sql << endl;
	    Query query5(conn);
	    query5.get_result(sql);
//			StoreQueryResult R3 = query5.store();
	    parentid = get_right_one(allids,query5);
	    //cout << "parentid2 " << parentid << endl;
	}
	route.push_back(parentid);
	/*
	 * add using the left right for the route
	 *
	 string left,right;
	 sql = "SELECT left_value,right_value FROM taxon WHERE taxon.taxon_id = "+nameid;
	 Query query5 = conn.query(sql);
	 StoreQueryResult R3 = query5.store();
	 left = R3[0][0].c_str();
	 right = R3[0][1].c_str();

	 vector<string> route;
	 sql = "SELECT taxon.taxon_id FROM taxon WHERE left_value < "+left;
	 sql += "AND right_value > ";
	 sql += right;
	 sql+= "AND left_value >= ";
	 sql+= cladeleft;
	 sql+= "AND right_value <=  ";
	 sql += claderight;
	 sql += "ORDER BY left_value DESC;";
	 Query query4 = conn.query(sql);
	 StoreQueryResult R2 = query4.store();
	 for(int j=0;j<R2.size();j++){
	 route.push_back(R2[j][0].c_str());
	 }
	 *
	 * end add
	 */
	vector<double> tdistance;
	for(int j=0;j<names.size();j++){
	    if(j!=i){
		sql = "SELECT ncbi_id FROM taxonomy WHERE name = '"+names[j]+"'  and name_class = 'scientific name';";
		Query query5(conn);
		query5.get_result(sql);
//				StoreQueryResult R3 = query5.store();
		string jnameid = get_right_one(allids,query5);
		double distance = 0;
		while((int)count(route.begin(),route.end(),jnameid) == 0){
		    distance += 1;
		    sql = "SELECT parent_ncbi_id FROM taxonomy WHERE ncbi_id = "+jnameid+" and name_class = 'scientific name';";
		    Query query6(conn);
		    query6.get_result(sql);
//					StoreQueryResult R4 = query6.store();
		    jnameid = get_right_one(allids,query6);
		}
		for(int k=0;k<route.size();k++){
		    if(route[k] == jnameid)
			distance += k;
		}
		tdistance.push_back(distance);
	    }else{
		tdistance.push_back(666666);
	    }
	}
	std::ostringstream stm;
	stm << i;
	numnames->insert( pair<string,string>(stm.str(),names[i]) );
	namesnum->insert( pair<string,string>(names[i],stm.str()) );
	numlist->push_back(tdistance);
	cout << "distances complete: "<< names[i] << endl;
    }
    //distances = list(tdistance)
}

//user tree one
void SQLiteProfiler::create_distances(vector<string> names,Tree * tree,map<string,string> * numnames
		,map<string,string> * namesnum, vector< vector<double> > * numlist){
    for(int i=0;i<names.size();i++){
	Node * nd1 = tree->getExternalNode(names[i]);
	vector<double> tdistance;
	for(int j=0;j<names.size();j++){
	    if (i == j){
		tdistance.push_back(666666);
	    }else{
		Node * nd2 = tree->getExternalNode(names[j]);
		double distance = get_distance_between_two_nodes(tree,nd1,nd2);
		tdistance.push_back(distance);
	    }
	}
	std::ostringstream stm;
	stm << i;
	numnames->insert( pair<string,string>(stm.str(),names[i]) );
	namesnum->insert( pair<string,string>(names[i],stm.str()) );
	numlist->push_back(tdistance);
    }
}

void SQLiteProfiler::get_shortest_distance_with_dicts(vector<string> names, map<string,string> numnames,map<string,string> namesnum,
		vector< vector<double> > numlist, string * shortestnameone, vector<string> * shortestnametwo){
    vector<double> distances;
    double shortestdistance = 10000;
    for(int i=0;i<names.size();i++){
	bool keepD = false;
	string nameid = namesnum[names[i]];
	vector<double> tdistance = numlist[atoi(nameid.c_str())];
	//for(int j=0;j<names.size();j++){
	//	cout << tdistance[j] << " " ;
	//}
	//cout << endl;
	for(int j=0;j<file_names.size();j++){
	    string nameid2 = namesnum[file_names[j]];
	    double distance = tdistance[atoi(nameid2.c_str())];
	    if(distance < shortestdistance && distance != 666666){
		shortestdistance = distance;
		*shortestnameone = names[i];
		keepD = true;
	    }
	}
	if (keepD == true){
	    distances = tdistance;
	}
    }
    shortestnametwo->clear();
    for(int j=0;j<file_names.size();j++){
	//int ct = (int) count(shortestnametwo->begin(),shortestnametwo->end(),file_names[j]);
	if(distances[j]==shortestdistance && distances[j]!=666666){
	    shortestnametwo->push_back(file_names[j]);
	    cout << "f " <<file_names[j]<<endl;
	}
    }
}

void SQLiteProfiler::clean_before_profile(string infile){
    string tfilen = infile;
    fix_bad_chars(tfilen);
    string cmd = "phyutility -clean 0.10 -in ";
    cmd += profilefoldername;
    cmd += tfilen;
    cmd += " -out ";
    cmd += profilefoldername;
    cmd += tfilen;

    cout << cmd << endl;
    FILE *fp = popen(cmd.c_str(), "r" );
    char buff[1000];
    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	string line(buff);
    }
    pclose( fp );
}

/*
 * TODO: fix the orphans
 */
void SQLiteProfiler::profile(map<string,string> numnames,map<string,string> namesnum,
		vector< vector<double> > numlist){
    cout << "writing everything to record.log" << endl;
    string recordname = profilefoldername+"/record.log";
    ofstream ofs(recordname.c_str());
    bool muscle = true;
    vector<string> profile_files;
    map<string,int> profile_files_dict;
    int profile_files_count = 0;
    vector<string> newnames(file_names.begin(),file_names.end());
    string last_profile_file; //for changing to FINAL.aln
    while (newnames.size() > 0){
	string shortestnameone = "";
	vector<string> * shortestnamestwo = new vector<string>();
	get_shortest_distance_with_dicts(newnames,numnames,namesnum,numlist, &shortestnameone, shortestnamestwo);
	cout << "shortestnameone " << shortestnameone << endl;
	for(int i=0;i<shortestnamestwo->size();i++)
	    cout << shortestnamestwo->at(i) << " ";
	cout << endl;
	vector<string>::iterator it;
	it = find (newnames.begin(), newnames.end(), shortestnameone);
	newnames.erase(it);
	//only one, simple case
	if(shortestnamestwo->size() == 1){
	    string firstfile = shortestnameone;
	    string secondfile;
	    bool already = false;
	    for(int i=0; i < profile_files.size(); i++){
		for(int j=0;j<shortestnamestwo->size();j++){
		    size_t found;
		    found=profile_files[i].find(shortestnamestwo->at(j));
		    if(found!=string::npos){
			already = true;
			secondfile = profile_files[i];
		    }
		}
	    }
	    if(already == false){
		secondfile = shortestnamestwo->at(0);
	    }
	    profile_files.push_back(firstfile+"__"+secondfile);
	    profile_files_dict.insert( pair<string,int>(firstfile+"__"+secondfile,profile_files_count) );
	    ofs << profile_files_count << "," << firstfile+"++"+secondfile << endl;
	    if(muscle == false){
		//if already == True:
		//	os.system("mafft-linsi --seed "+directory+str(profile_files_dict[secondfile])+" "+directory+firstfile+" > "+directory+str(profile_files_dict[firstfile+"__"+secondfile]))
		//else:
		//	os.system("mafft-linsi --seed "+directory+secondfile+" "+directory+firstfile+" > "+directory+str(profile_files_dict[firstfile+"__"+secondfile]))
	    }else{
		if(already == true){
		    string cmd = "muscle -profile -in1 ";
		    cmd += profilefoldername;
		    std::string ts;
		    std::stringstream tout;
		    tout << profile_files_dict[secondfile];
		    cout << secondfile << " " << profile_files_dict[secondfile] << endl;
		    ts = tout.str();
		    cmd += ts;
		    cmd += " -in2 ";
		    cmd += profilefoldername;
		    string tfilen = firstfile;
		    fix_bad_chars(tfilen);
		    cmd += tfilen;
		    //cmd += firstfile;//fix space
		    cmd += " -out ";
		    cmd += profilefoldername;
		    std::string ts2;
		    std::stringstream tout2;
		    tout2 << profile_files_dict[firstfile+"__"+secondfile];
		    ts2 = tout2.str();
		    cmd += ts2;
		    cout << secondfile << " " << profile_files_dict[firstfile+"__"+secondfile] << endl;
		    cout << "aligning1" << endl;
		    cout << cmd << endl;
		    FILE *fp = popen(cmd.c_str(), "r" );
		    char buff[1000];
		    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
			string line(buff);
		    }
		    pclose( fp );
		    //clean
		    clean_before_profile(ts2);
		    last_profile_file = ts2;
		}else{
		    string cmd = "muscle -profile -in1 ";
		    cmd += profilefoldername;
		    string tfilen2 = secondfile;
		    fix_bad_chars(tfilen2);
		    cmd += tfilen2;
		    //cmd += secondfile;//fix space
		    cmd += " -in2 ";
		    cmd += profilefoldername;
		    string tfilen = firstfile;
		    fix_bad_chars(tfilen);
		    cmd += tfilen;
		    //cmd += firstfile;//fix space here
		    cmd += " -out ";
		    cmd += profilefoldername;
		    std::string s;
		    std::stringstream out;
		    out << profile_files_dict[firstfile+"__"+secondfile];
		    s = out.str();
		    cmd += s;
		    cout << "aligning2" << endl;
		    cout <<cmd << endl;
		    FILE *fp = popen(cmd.c_str(), "r" );
		    char buff[1000];
		    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
			string line(buff);
		    }
		    pclose( fp );
		    //clean
		    clean_before_profile(s);
		    last_profile_file = s;
		}
	    }
	    cout << "moving on " << profile_files_count << endl;
	    profile_files_count += 1;
	    if(already == false){
		vector<string>::iterator it;
		it = find (newnames.begin(), newnames.end(), secondfile);
		newnames.erase(it);
	    }else{
		vector<string>::iterator it;
		it = find (profile_files.begin(), profile_files.end(), secondfile);
		//cout <<profile_files.size()<<endl;
		profile_files.erase(it);
		//cout <<profile_files.size()<<endl;
	    }
	}else{
	    string firstfile = shortestnameone;
	    string secondfile;
	    bool already = false;
	    for(int i=0; i < profile_files.size(); i++){
		for(int j=0;j<shortestnamestwo->size();j++){
		    size_t found;
		    found=profile_files[i].find(shortestnamestwo->at(j));
		    if(found!=string::npos){
			already = true;
			secondfile = profile_files[i];
		    }
		}
	    }
	    if (already == false){
		secondfile = shortestnamestwo->at(0);
		string bestsn = shortestnamestwo->at(0);
		double bestscore = 0;
		for(int i=0;i<shortestnamestwo->size();i++){
		    if(muscle == false){
			//os.system("mafft-linsi --seed "+directory+sn+" "+directory+firstfile+" > "+directory+firstfile+"__"+sn)
		    }else{
			string cmd = "muscle -profile -in1 ";
			cmd += profilefoldername;
			string tfilen2(shortestnamestwo->at(i));
			fix_bad_chars(tfilen2);
			cmd += tfilen2;
			//cmd += shortestnamestwo->at(i);//fix space here
			cmd += " -in2 ";
			cmd += profilefoldername;
			string tfilen = firstfile;
			fix_bad_chars(tfilen);
			cmd += tfilen;
			//cmd += firstfile;//fix space here
			cmd += " -out ";
			cmd += profilefoldername;
			cmd += tfilen+"__"+tfilen2;
			cout << "aligning many1" << endl;
			cout << cmd << endl;
			FILE *fp = popen(cmd.c_str(), "r" );
			char buff[1000];
			while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
			    string line(buff);
			}
			pclose( fp );
		    }
		    profile_files_count += 1;
		    string cmd = "muscle -spscore ";
		    cmd += profilefoldername;
		    string tfilen = firstfile;
		    fix_bad_chars(tfilen);
		    string tfilen2 = shortestnamestwo->at(i);
		    fix_bad_chars(tfilen2);
		    cmd += tfilen+"__"+tfilen2;
		    cmd += " -log prlog";
		    cout << "aligning many2" << endl;
		    cout << cmd << endl;
		    FILE *fp = popen(cmd.c_str(), "r" );
		    char buff[1000];
		    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
			string line(buff);
		    }
		    pclose( fp );
		    //read file
		    double score = 0;
		    ifstream ifs("prlog");
		    string line;
		    while(getline(ifs,line)){
			TrimSpaces(line);
			size_t se = line.find("="); // Find the first character position from reverse af
			// if all spaces or empty return an empty string
			if (se != string::npos){
			    vector<string> tokens;
			    string del("SP=");
			    Tokenize(line, tokens, del);
			    score = atof(tokens[2].c_str());
			}
		    }
		    ifs.close();
		    if(score > bestscore){
			bestscore = score;
			bestsn = shortestnamestwo->at(i);
		    }
		}
		secondfile = bestsn;
		profile_files.push_back(firstfile+"__"+secondfile);
		profile_files_dict.insert( pair<string,int>(firstfile+"__"+secondfile,profile_files_count) );
		ofs << profile_files_count << "," << firstfile+"++"+secondfile << endl;
		string cmd = "cp ";
		cmd += profilefoldername;
		string tfilen = firstfile;
		fix_bad_chars(tfilen);
		string tfilen2 = secondfile;
		fix_bad_chars(tfilen2);
		cmd += tfilen+"__"+tfilen2;
		cmd += " ";
		cmd += profilefoldername;
		std::string s;
		std::stringstream out;
		out << profile_files_dict[firstfile+"__"+secondfile];
		s = out.str();
		cmd += s;
		cout << "cp many3" << endl;
		cout << cmd << endl;
		FILE *fp = popen(cmd.c_str(), "r" );
		char buff[1000];
		while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
		    string line(buff);
		}
		pclose( fp );
		//clean
		clean_before_profile(s);
		last_profile_file = s;
		profile_files_count += 1;
		vector<string>::iterator it;
		it = find (newnames.begin(), newnames.end(), secondfile);
		newnames.erase(it);
	    }else{
		//mafft profile
		vector<string>::iterator it;
		it = find (profile_files.begin(), profile_files.end(), secondfile);
		profile_files.erase(it);
		profile_files.push_back(firstfile+"__"+secondfile);
		profile_files_dict.insert( pair<string,int>(firstfile+"__"+secondfile,profile_files_count) );
		ofs << profile_files_count << "," << firstfile+"++"+secondfile << endl;
		if (muscle == false){
		    //	os.system("mafft-linsi --seed "+directory+str(profile_files_dict[secondfile])+" "+directory+firstfile+" > "+directory+str(profile_files_dict[firstfile+"__"+secondfile]))
		}else{
		    //	os.system("muscle -profile -in1 "+directory+str(profile_files_dict[secondfile])+" -in2 "+directory+firstfile+" -out "+directory+str(profile_files_dict[firstfile+"__"+secondfile]))
		    string cmd = "muscle -profile -in1 ";
		    cmd += profilefoldername;
		    std::string s;
		    std::stringstream out;
		    out << profile_files_dict[secondfile];
		    s = out.str();
		    cmd += s;
		    cmd += " -in2 ";
		    cmd += profilefoldername;
		    string tfilen = firstfile;
		    fix_bad_chars(tfilen);
		    cmd += tfilen;
		    //cmd += firstfile; //fix space here
		    cmd += " -out ";
		    cmd += profilefoldername;
		    std::string s2;
		    std::stringstream out2;
		    out2 << profile_files_dict[firstfile+"__"+secondfile];
		    s2 = out2.str();
		    cmd += s2;
		    cout << "aligning3" << endl;
		    cout << cmd << endl;
		    FILE *fp = popen(cmd.c_str(), "r" );
		    char buff[1000];
		    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
			string line(buff);
		    }
		    pclose( fp );
		    //clean
		    clean_before_profile(s2);
		    last_profile_file = s2;
		}
		profile_files_count += 1;
	    }
	}
	delete shortestnamestwo;
    }
    for(int i=0; i < profile_files.size(); i++){
	cout << profile_files[i] << " ++ "<< profile_files_dict[profile_files[i]] << endl;
    }
    while(profile_files.size() > 1){
	string shortestnameone = "";
	vector<string> * shortestnamestwo = new vector<string>();
	get_shortest_distance_with_dicts(file_names,numnames,namesnum,numlist, &shortestnameone, shortestnamestwo);
	vector<string> tokens;
	string del("_");
	Tokenize(profile_files[0], tokens, del);
	string test1 = tokens[0].c_str();
	string numtest1 = namesnum[test1];
	string bestfile = "";
	double bestsc = 1000000;
	for(int i=1;i<profile_files.size();i++){
	    double tdis = 100000;
	    tokens.clear();
	    Tokenize(profile_files[i], tokens, del);
	    for(int j=0;j<tokens.size();j++){
		//cout << tokens[j] << endl;
		//cout << numtest1 << "==" << namesnum[tokens[j]] <<endl;
		//cout << numlist[atoi(numtest1.c_str())][atoi(namesnum[tokens[j]].c_str())] << endl;
		if (tdis > numlist[atoi(numtest1.c_str())][atoi(namesnum[tokens[j]].c_str())]){
		    tdis = numlist[atoi(numtest1.c_str())][atoi(namesnum[tokens[j]].c_str())];
		}
	    }
	    if(tdis < bestsc){
		bestsc = tdis;
		bestfile = profile_files[i];
	    }
	}
	string firstfile = profile_files[0];
	string secondfile = bestfile;
	vector<string>::iterator it;
	it = find (profile_files.begin(), profile_files.end(), firstfile);
	profile_files.erase(it);
	it = find (profile_files.begin(), profile_files.end(), secondfile);
	profile_files.erase(it);

	profile_files.push_back(firstfile+"__"+secondfile);
	profile_files_dict.insert( pair<string,int>(firstfile+"__"+secondfile,profile_files_count) );
	ofs << profile_files_count << "," << firstfile+"++"+secondfile <<endl;
	if (muscle == false){
	    //	os.system("mafft-linsi --seed "+directory+str(profile_files_dict[secondfile])+" "+directory+firstfile+" > "+directory+str(profile_files_dict[firstfile+"__"+secondfile]))
	}else{
	    //	os.system("muscle -profile -in1 "+directory+str(profile_files_dict[secondfile])+" -in2 "+directory+firstfile+" -out "+directory+str(profile_files_dict[firstfile+"__"+secondfile]))
	    string cmd = "muscle -profile -in1 ";
	    cmd += profilefoldername;
	    std::string s1;
	    std::stringstream out1;
	    out1 << profile_files_dict[firstfile];
	    s1 = out1.str();
	    cmd += s1;
	    cmd += " -in2 ";
	    cmd += profilefoldername;
	    std::string s3;
	    std::stringstream out3;
	    out3 << profile_files_dict[secondfile];
	    s3 = out3.str();
	    cmd += s3;
	    cmd += " -out ";
	    cmd += profilefoldername;
	    std::string s2;
	    std::stringstream out2;
	    out2 << profile_files_dict[firstfile+"__"+secondfile];
	    s2 = out2.str();
	    cmd += s2;
	    cout << "aligningfinal" << endl;
	    cout << cmd << endl;
	    FILE *fp = popen(cmd.c_str(), "r" );
	    char buff[1000];
	    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
		string line(buff);
	    }
	    pclose( fp );
	    //clean
	    clean_before_profile(s2);
	    last_profile_file = s2;
	}
	profile_files_count += 1;
    }

    //do the orphans
    //print profile_files
    //print profile_files_dict[profile_files[0]]
    /*
     * fix orphans
     */
    //if useorphan == True:
    //	os.system("mafft-linsi --seed "+orphanfile+" "+directory+str(profile_files_dict[profile_files[0]])+" > "+directory+"FINAL.aln")
    //else:
    //	os.system("cp "+directory+profile_files[0]+" "+directory+"FINAL.aln")
    copy_final_file(last_profile_file);
    ofs.close();
}

/*
 * The idea here is to run through the profilerun (and therefore profilekey)
 * rerunning all the flaggedprofiles (in the vector)
 *
 * They are already premapped so you are just rerunning the ones that have updated seqs
 */
void SQLiteProfiler::update_profile(){
	for(unsigned int i=0;i<profilerun.size();i++){
		cout << count(flaggedprofiles.begin(),flaggedprofiles.end(),profilerun[i])  << " " << profilerun[i] << endl;
		if(count(flaggedprofiles.begin(),flaggedprofiles.end(),profilerun[i]) > 0){
			cout << "reruning " << profilerun[i] << endl;
			string profile1 = get_profilekey_value(profilekey[profilerun[i]][0]);
			string profile2 = get_profilekey_value(profilekey[profilerun[i]][1]);
			cout << profile1 << " " << profile2 <<endl;
			string cmd = "muscle -profile -in1 ";
			cmd += profilefoldername;
			std::string ts;
			std::stringstream tout;
			tout << profile1;
			ts = tout.str();
			cmd += ts;
			cmd += " -in2 ";
			cmd += profilefoldername;
			string tfilen = profile2;
			fix_bad_chars(tfilen);
			cmd += tfilen;
			cmd += " -out ";
			cmd += profilefoldername;
			std::string ts2;
			std::stringstream tout2;
			tout2 << profilerun[i];
			ts2 = tout2.str();
			cmd += ts2;
			cout << "aligning" << endl;
			cout << cmd << endl;
			FILE *fp = popen(cmd.c_str(), "r" );
			char buff[1000];
			while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
				string line(buff);
			}
			pclose( fp );
			//clean
			clean_before_profile(ts2);
		}else{
			cout << "skipping " << profilerun[i] << endl;
		}
	}
	if(profilerun.size() > 0){
		copy_final_file(profilerun[profilerun.size()-1]);
	}else{
		copy_final_file(updatedfiles[0]); //should mean that there is only one file
	}
}

/*
 * copy final file to FINAL.aln
 */

void SQLiteProfiler::copy_final_file(string filename){
    string tfilen = filename;
    fix_bad_chars(tfilen);
    string cmd = "cp ";
    cmd += profilefoldername;
    cmd += tfilen;
    cmd += " ";
    cmd += profilefoldername;
    cmd += "FINAL.aln";
    cout << cmd << endl;
    FILE *fp = popen(cmd.c_str(), "r" );
    char buff[1000];
    while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	string line(buff);
    }
    pclose( fp );
}

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
	    /*
	     * record the id as the first token
	     * and the numbers below as the numbers
	     */
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
    /*
     * calculate the means
     */
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
    /*
     * calculate the means
     */
    vector<double> mns(numseqs);
    for(int i=0;i<numseqs;i++){
	//TODO : there was a pointer problem here
	mns[i] = allmeans[sequences->at(i).get_id()];
    }
    double sd = stdev(mns);
    double mn = mean(mns);
    double dev = sd*1.5+mn;//obviously changeable
    /*
     * read in the fasta file
     */
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
void SQLiteProfiler::rename_final_alignment(string which){
    /*
     * read in the final file
     */
    FastaUtil seqreader;
    vector<Sequence> * sequences = new vector<Sequence>();
    seqreader.readFile(profilefoldername+which, *sequences);
    /*
     * read in the table
     */
    string line;
    string tblname = gene_name+".gi";
    ifstream pfile (tblname.c_str());
    vector<string> tokens;
    bool first = true;
    map<string,string> db_to_ncbi;
    map<string,string> db_to_name;
    if (pfile.is_open()){
	while (! pfile.eof() ){
	    if (first == false){
		getline (pfile,line);
		string del("\t");
		tokens.clear();
		Tokenize(line, tokens, del);
		if(tokens.size() >= 1){
		    db_to_ncbi[tokens[0]] = tokens[1];
		    string astr = get_name_from_tax_id(tokens[0]);
		    fix_bad_chars_for_seq_names(astr);
		    db_to_name[tokens[0]] = astr;
		}
	    }else{
		getline (pfile,line);
		first = false;
	    }
	}
    }
    /*
     * write the final file
     */
    vector<Sequence> sc1;
    FastaUtil seqwriter;
    string fn1 = profilefoldername+which+".rn";
    for(int i=0;i<sequences->size();i++){
	//TODO : there was a pointer problem here
	string tname = db_to_name[sequences->at(i).get_id()];
	Sequence ts = sequences->at(i);
	ts.set_id(tname);
	sc1.push_back(ts);
    }
    seqwriter.writeFileFromVector(fn1,sc1);
    delete sequences;
}


#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <stack>
#include <algorithm>

#include "genedb.h"
#include "libsqlitewrapped.h"
#include "utils.h"
#include "fasta_util.h"
#include "sequence.h"

using namespace std;

template <class T>
inline std::string to_string (const T& t){
    std::stringstream ss;
    ss << t;
    return ss.str();
}

GeneDB::GeneDB(){}

GeneDB::GeneDB(string nm):name(nm){}

/*
 * Should check to see if the database exists, if not 
 * then it will initialize it. If so then it will make
 * sure that the structure is correct.
 */
/*
 * SQL structure as of 0.3a
 * 
 * alignments
 * -
 * 
 * sequences
 * -
 *
 * sequence_alignment_map
 * -
 * 
 *
 * alignment_alignment_map
 * -
 * 
 * record
 * -
 */
void GeneDB::initialize(bool overwrite){
    //check to see if the file exists
    ifstream ifile(name.c_str());
    if (ifile && overwrite == false){
	cerr << name << " exists" <<endl;
	//TODO: add a check to make sure that it has the right structure
    }else{
	if(ifile && overwrite == true){//initializing the database with the SQL structure
	    remove(name.c_str());
	}
	Database conn(name);
	cout << "connected to " << name << endl;
	Query query(conn);
	query.execute("create table alignments(id INTEGER PRIMARY KEY, alignname VARCHAR(255),updated INTEGER);");
	query.execute("create index alignments_alignname on alignments(alignname);");
	query.execute("create table profile_alignments(id INTEGER PRIMARY KEY, child1 INTEGER, child2 INTEGER, name VARCHAR(255));");
//make index
	query.execute("create table sequence_profile_map(id INTEGER PRIMARY KEY, sequence_id INTEGER, profile_id INTEGER, sequence TEXT);");
//make index
	query.execute("create table sequences(id INTEGER PRIMARY KEY, ncbi_id INTEGER, accession VARCHAR(255), name VARCHAR(255), sequence TEXT);");
	query.execute("create index sequences_ncbi_id on sequences(ncbi_id);");
	query.execute("create index sequences_name on sequences(name);");
	query.execute("create table sequence_alignment_map(id INTEGER PRIMARY KEY, sequence_id INTEGER, alignment_id INTEGER, sequence TEXT);");
	query.execute("create index sequence_alignment_map_sequence_id sequence_alignment_map(sequence_id);");
	query.execute("create index sequence_alignment_map_alignment_id sequence_alignment_map(alignment_id);");
	cout << "database created"<<endl;
    }
}

void GeneDB::add_user_seqs_to_db(vector<Sequence> * user_seqs){
    cout << "adding user sequences to database: " << name << endl;
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    for(int i=0;i<user_seqs->size();i++){
	string sql = "insert into sequences (ncbi_id,accession,name,sequence) values (";
	sql += user_seqs->at(i).get_ncbi_tax_id()+",0,'";
	sql += user_seqs->at(i).get_id()+"','";
	sql += user_seqs->at(i).get_sequence()+"');";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
	int sid = sqlite3_last_insert_rowid(conn);
	user_seqs->at(i).set_sqlite_id(sid);
    }
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);
}

void GeneDB::add_seqs_to_db(vector<Sequence> * keep_seqs){
    cout << "adding sequences to database: " << name << endl;
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    for(int i=0;i<keep_seqs->size();i++){
	string sql = "insert into sequences (ncbi_id,accession,name,sequence) values (";
	sql += keep_seqs->at(i).get_ncbi_tax_id()+",'";
	sql += keep_seqs->at(i).get_ncbi_gi_id()+"','";
	sql += keep_seqs->at(i).get_name()+"','";
	sql += keep_seqs->at(i).get_sequence()+"');";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
	int sid = sqlite3_last_insert_rowid(conn);
	keep_seqs->at(i).set_sqlite_id(sid);
    }
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);
}

//TODO: add the reading of the actual alignment file
int GeneDB::add_alignment(string filename, vector<Sequence> * dbseqs, vector<Sequence> * userseqs){
    cout << "adding alignment to database: " << name << endl;
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string sql = "insert into alignments (alignname,updated) values ('";
    sql += filename+"',0);";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    int alignid = sqlite3_last_insert_rowid(conn);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    string alignids = to_string(alignid);
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    for(int i=0;i<dbseqs->size();i++){
	string sql = "insert into sequence_alignment_map (sequence_id, alignment_id,sequence) values (";
	string dbsid = to_string(dbseqs->at(i).get_sqlite_id());
	sql += dbsid+",";
	sql += alignids+",'";
	sql += dbseqs->at(i).get_aligned_seq()+"');";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    }
    for(int i=0;i<userseqs->size();i++){
	string sql = "insert into sequence_alignment_map (sequence_id, alignment_id,sequence) values (";
	string dbsid = to_string(userseqs->at(i).get_sqlite_id());
	sql += dbsid+",";
	sql += alignids+",'";
	sql += userseqs->at(i).get_aligned_seq()+"');";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    }
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);
    return alignid;
}

void GeneDB::toggle_alignment_update(int alignid){
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string sql = "update alignments set updated = 1";
    sql += " where id = ";
    sql += to_string(alignid)+";";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
}

void GeneDB::remove_alignment(int alignid){
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string sql = "delete from sequence_alignment_map where alignment_id =";
    string alignids = to_string(alignid);
    sql += alignids+";";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);

    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    sql = "delete from alignments where id =";
    sql += alignids+";";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
}

void GeneDB::remove_profile_alignments(){
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string sql = "delete from sequence_profile_map;";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);

    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    sql = "delete from profile_alignments;";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
}

/*
 * these would be the childless alignments, the first ones
 */
void GeneDB::get_first_profile_alignments(vector<string> & names){
    string sql = "select name from profile_alignments where child1 = 0 and child2 = 0;";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    while(query.fetch_row()){
	names.push_back(to_string(query.getstr()));
    }
    query.free_result();
}

int GeneDB::get_alignment_id_by_name(string alignname){
    string sql = "select id from alignments where alignname = '";
    sql += alignname+"';";

    Database conn(name);
    Query query(conn);
    int alignid;
    query.get_result(sql);
    while(query.fetch_row()){
	alignid = atoi(to_string(query.getval()).c_str());
    }
    query.free_result();
    return alignid;
}

void GeneDB::remove_alignment_by_name(string alignname){
    int alignid = get_alignment_id_by_name(alignname);
    remove_alignment(alignid);
}

/*
 * uses the get_sqlite_id for the sqlite id
 */
void GeneDB::add_seq_to_alignment(int alignid, Sequence inseq){
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string alignids = to_string(alignid);
    string sql = "insert into sequence_alignment_map (sequence_id, alignment_id,sequence) values (";
    string dbsid = to_string(inseq.get_sqlite_id());
    sql += dbsid+",";
    sql += alignids+",'";
    sql += inseq.get_aligned_seq()+"');";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
}

/*
 * seqs is returned with the id being the sqlite id and the seq being the 
 * current aligned seq
 */
void GeneDB::get_align_seqs(int alignid, vector<Sequence> & seqs){
    
    string alignids = to_string(alignid);
    string sql = "select sequence_id,sequence from sequence_alignment_map where alignment_id = ";
    sql += alignids+";";

    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    while(query.fetch_row()){
	string id;
	string seq;
	id = to_string(query.getval());
	seq = to_string(query.getval());
	Sequence tseq(id,seq);
	seqs.push_back(tseq);
    }
    query.free_result();
}

/*
 * seqs is returned with the id being the sqlite id and the seq being the 
 * original unaligned seq
 */
void GeneDB::get_align_seqs_unaligned(int alignid, vector<Sequence> & seqs){
    string alignids = to_string(alignid);
    string sql = "select sequence_alignment_map.sequence_id,sequences.sequence from sequence_alignment_map,sequences where sequence_alignment_map.sequence_id=sequences.id and sequence_alignment_map.alignment_id = ";
    sql += alignids+";";

    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    while(query.fetch_row()){
	string id;
	string seq;
	id = to_string(query.getval());
	seq = to_string(query.getstr());
	Sequence tseq(id,seq);
	seqs.push_back(tseq);
    }
    query.free_result();
}

/*
 * seqs is returned with the id being the ncbi_id and the seq being the 
 * original unaligned seq and the other elements correctly defined
 * the alignname is the name in the filename column
 */
void GeneDB::get_align_seq_unaligned_fully_initialized(string alignname,vector<Sequence> & seqs){
    int alignid = get_alignment_id_by_name(alignname);
    string sql = "select sequence_alignment_map.sequence_id,sequences.ncbi_id,sequences.accession,sequences.name,sequences.sequence from sequence_alignment_map,sequences where sequence_alignment_map.sequence_id=sequences.id and sequence_alignment_map.alignment_id = ";
    sql += to_string(alignid)+";";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    while(query.fetch_row()){
	string id;
	string ncbi_id;
	string acc;
	string name;
	string seq;
	id = to_string(query.getval());
	ncbi_id = to_string(query.getval());
	acc = to_string(query.getstr());
	name = to_string(query.getstr());
	seq = to_string(query.getstr());
	Sequence tseq(ncbi_id,seq);
	if (name.find("user_") != string::npos){//user id is username if user seq
	    tseq.set_id(name);
	}
	tseq.set_ncbi_tax_id(ncbi_id);
	tseq.set_ncbi_gi_id(acc);
	tseq.set_name(name);
	tseq.set_sqlite_id(atoi(id.c_str()));
	seqs.push_back(tseq);
    }
    query.free_result();
}

/*
 * uses the get_sqlite_id and the get_aligned_seq
 */
void GeneDB::update_align_seqs(int alignid,vector<Sequence> & seqs){
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    for(int i=0;i<seqs.size();i++){
	string sql = "update sequence_alignment_map set sequence = '";
	sql += seqs[i].get_aligned_seq()+"'";
	sql += "' where alignment_id = ";
	sql += to_string(alignid)+" and sequence_id = ";
	sql += to_string(seqs[i].get_sqlite_id())+";";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    }
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
}

/*
 * ncbi id will be the id and the comment
 * name will be the name
 * seq will be the original seq
 * sqlite_id will be the sqlite_id
 */
void GeneDB::get_all_sequences(vector<Sequence> & seqs){
    string sql = "select * from sequences;";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    while(query.fetch_row()){
	string id;
	string ncbi_id;
	string acc;
	string name;
	string seq;
	id = to_string(query.getval());
	ncbi_id = to_string(query.getval());
	acc = to_string(query.getstr());
	name = to_string(query.getstr());
	seq = to_string(query.getstr());
	Sequence tseq(ncbi_id,seq);
	if (name.find("user_") != string::npos){//user id is username if user seq
	    tseq.set_id(name);
	}
	tseq.set_ncbi_tax_id(ncbi_id);
	tseq.set_ncbi_gi_id(acc);
	tseq.set_sqlite_id(atoi(id.c_str()));
	tseq.set_name(name);
	seqs.push_back(tseq);
    }
    query.free_result();
}

void GeneDB::get_alignment_names(vector<string> & names){
    string sql = "select alignname from alignments;";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    while(query.fetch_row()){
	string name;
	name = to_string(query.getstr());
	names.push_back(name);
    }
    query.free_result();
}

void GeneDB::get_alignment_nums(vector<int> & nums){
    string sql = "select id from alignments;";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    while(query.fetch_row()){
	int num;
	num = query.getval();
	nums.push_back(num);
    }
    query.free_result();
}

/* 
 * typically used for updating and should be called after moving over updated 
 */
void GeneDB::get_profile_alignment_nums(vector<int> & nums){
    string sql = "select id from profile_alignments where child1 = 0 and child2 = 0;";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    while(query.fetch_row()){
	int num;
	num = query.getval();
	nums.push_back(num);
    }
    query.free_result();
}

void GeneDB::copy_alignments_to_first_profiles(map<int, string> & profile_id_name_map){
    string sql = "select id,alignname from alignments;";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    vector<int> aid;
    vector<string> names;
    while(query.fetch_row()){
	int tid;
	string name;
	tid = query.getval();
	name = to_string(query.getstr());
	aid.push_back(tid);
	names.push_back(name);
    }
    query.free_result();
    //add the new alignment
    
    
    for (int i=0;i<names.size();i++){
	sqlite3 *conn2;
	int rc = sqlite3_open(name.c_str(), &conn2);
	char *zErrMsg = 0;
	sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
	string sql = "insert into profile_alignments (child1,child2,name) values (0,0,'";
	sql += names[i]+"');";
	rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
	int pid = sqlite3_last_insert_rowid(conn2);
	sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
	profile_id_name_map[pid] = names[i];

	sql = "select sequence_id,sequence from sequence_alignment_map where alignment_id = ";
	sql += to_string(aid[i])+";";
	Query query2(conn);
	query2.get_result(sql);
	vector<int> sid;
	vector<string> seqs;
	while(query2.fetch_row()){
	    sid.push_back(query2.getval());
	    seqs.push_back(to_string(query2.getstr()));
	}
	query2.free_result();

	//insert the seqs into profile_alignment sequences
	sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
	for(int j=0;j<sid.size();j++){
	    sql = "insert into sequence_profile_map (sequence_id,profile_id,sequence) values (";
	    sql += to_string(sid[j])+",";
	    sql += to_string(pid)+",'";
	    sql += seqs[j]+"');";
	    rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
	}
	sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
	sqlite3_close(conn2);

    }
}

void GeneDB::copy_alignments_to_first_profiles_updated(map<int, string> & profile_id_name_map,vector<int>& updatedprofsnums){
    string sql = "select id,alignname from alignments where updated = 1;";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    vector<int> aid;
    vector<string> names;
    while(query.fetch_row()){
	int iid;
	string name;
	iid = query.getval();
	name = to_string(query.getstr());
	aid.push_back(iid);
	names.push_back(name);
    }
    query.free_result();
    //add the new alignment
    
    
    for (int i=0;i<names.size();i++){
	sqlite3 *conn2;
	int rc = sqlite3_open(name.c_str(), &conn2);
	char *zErrMsg = 0;
	sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
	string sql = "insert into profile_alignments (child1,child2,name) values (0,0,'";
	sql += names[i]+"');";
	rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
	int pid = sqlite3_last_insert_rowid(conn2);
	sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
	updatedprofsnums.push_back(pid);

	sql = "select sequence_id,sequence from sequence_alignment_map where alignment_id = ";
	sql += to_string(aid[i])+";";
	Query query2(conn);
	query2.get_result(sql);
	vector<int> sid;
	vector<string> seqs;
	while(query2.fetch_row()){
	    sid.push_back(query2.getval());
	    seqs.push_back(to_string(query2.getstr()));
	}
	query2.free_result();

	//insert the seqs into profile_alignment sequences
	sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
	for(int j=0;j<sid.size();j++){
	    sql = "insert into sequence_profile_map (sequence_id,profile_id,sequence) values (";
	    sql += to_string(sid[j])+",";
	    sql += to_string(pid)+",'";
	    sql += seqs[j]+"');";
	    rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
	}
	sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
	sqlite3_close(conn2);
    }
    sql = "select id,name from profile_alignments where name != 0;";
    Query query3(conn);
    query3.get_result(sql);
    while(query3.fetch_row()){
	int iid;
	string name;
	iid = query3.getval();
	name = to_string(query3.getstr());
	profile_id_name_map[iid] = name;
    }
    query3.free_result();

}

int GeneDB::get_deepest_profile_for_alignment(int alignid){
    //get the deepest alignment in the profile alignment table that includes the alignid
    stack<int> retstack;
    retstack.push(alignid);
    int finalret = -1;
    while(!retstack.empty()){
	int curid = retstack.top();
	retstack.pop();
	string sql = "select id from profile_alignments where child1 = ";
	sql += to_string(curid) +" or child2 = ";
	sql += to_string(curid) +";";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while(query.fetch_row()){
	    int tid;
	    string name;
	    tid = query.getval();
	    if (tid > finalret){
		finalret = tid;
	    } 
	    retstack.push(tid);
	}
	query.free_result();
    }
    return finalret;
}

int GeneDB::add_profile_alignment(int child_id1, int child_id2){
    sqlite3 *conn2;
    int rc = sqlite3_open(name.c_str(), &conn2);
    char *zErrMsg = 0;
    sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string sql = "insert into profile_alignments (child1,child2,name) values (";
    sql += to_string(child_id1)+",";
    sql += to_string(child_id2)+",'";
    sql += to_string(0)+"');";//not an ncbi 
    rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
//    cout << sql << " " << rc << endl;
    int pid = sqlite3_last_insert_rowid(conn2);
    sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
    return pid;
}

void GeneDB::write_profile_alignment_to_file(int alignid, string filename){
    string sql = "select sequence_id,sequence from sequence_profile_map where profile_id = ";
    sql += to_string(alignid) +";";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    vector<Sequence> seqs;
    while(query.fetch_row()){
	string name;
	string seq;
	name = to_string(query.getval());
	seq = to_string(query.getstr());
	Sequence tseq (name,seq);
	seqs.push_back(tseq);
    }
    query.free_result();
    FastaUtil fu;
    fu.writeFileFromVector(filename,seqs);
}

void GeneDB::write_profile_alignment_with_names_to_file(int alignid, string filename,bool ncbi){
    string sql = "select sequence_profile_map.sequence_id,sequences.ncbi_id,sequences.name,sequence_profile_map.sequence from sequence_profile_map,sequences where sequence_profile_map.sequence_id=sequences.id and sequence_profile_map.profile_id = ";
    sql += to_string(alignid) +";";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    vector<Sequence> seqs;
    while(query.fetch_row()){
	string sid;
	string ncbiid;
	string name;
	string seq;
	sid = to_string(query.getval());
	ncbiid = to_string(query.getval());
	name = to_string(query.getstr());
	seq = to_string(query.getstr());
	if (ncbi && ncbiid != "0")
	    name = ncbiid;
	Sequence tseq (name,seq);
	seqs.push_back(tseq);
    }
    query.free_result();
    FastaUtil fu;
    fu.writeFileFromVector(filename,seqs);
}


/*
 * assumes that the seq id is the sqlite id and the seq is the aligned seq
 */
void GeneDB::add_sequences_for_profile_alignment(int profilealignid,vector<Sequence> & seqs){
    sqlite3 *conn2;
    int rc = sqlite3_open(name.c_str(), &conn2);
    char *zErrMsg = 0;
    sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
    for(int i=0;i<seqs.size();i++){
	string sql = "insert into sequence_profile_map (sequence_id,profile_id,sequence) values (";
	sql += to_string(seqs[i].get_id())+",";
	sql += to_string(profilealignid)+",'";
	sql += to_string(seqs[i].get_sequence())+"');";//not an ncbi 
	rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
    }
    sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
}


void GeneDB::get_updated_profs_names_delete_old(vector<string> & updatedprofs,vector<int> & updatedprofsnums, vector<int>&notupdatedprofnums){
    //get the updated names
    string sql = "select alignments.id,alignments.alignname,profile_alignments.id from alignments,profile_alignments where alignments.alignname=profile_alignments.name and alignments.updated = 1;";
    Database conn(name);
    Query query(conn);
    query.get_result(sql);
    stack<int> retstack;
    vector<int> todelete;
    while(query.fetch_row()){
	int sid;
	string name;
	int aid;
	sid = query.getval();
	name = to_string(query.getstr());
	aid = query.getval();
	updatedprofs.push_back(name);
	updatedprofsnums.push_back(aid);
	retstack.push(aid);
	todelete.push_back(int(aid));
    }
    query.free_result();
    sql = "select id,name from profile_alignments where child1 = 0 and child2 = 0;";
    Query queryt(conn);
    queryt.get_result(sql);
    while(queryt.fetch_row()){
	int sid;
	string name;
	sid = queryt.getval();
	name = to_string(queryt.getstr());
//	cout << "notupdatedsearch: " << sid << " " << name << endl;
	if(count(updatedprofs.begin(),updatedprofs.end(),name)==0){
	    notupdatedprofnums.push_back(int(sid));
//	    cout << "sid: " << sid << endl;
	}
    }
    queryt.free_result();

    //get all the profiles involving these names
    int finalret = -1;
    while(!retstack.empty()){
	int curid = retstack.top();
	retstack.pop();
	string sql = "select id from profile_alignments where child1 = ";
	sql += to_string(curid) +" or child2 = ";
	sql += to_string(curid) +";";
	Database conn(name);
	Query query(conn);
	query.get_result(sql);
	while(query.fetch_row()){
	    int tid;
	    string name;
	    tid = query.getval();
	    if (count(todelete.begin(),todelete.end(),int(tid)) ==0)
		todelete.push_back(tid);
	    retstack.push(tid);
	}
	query.free_result();
    }

    //remove all the profiles that involve the names
    //TODO: could get the process
    for(int i=0;i<todelete.size();i++){
	sqlite3 *conn2;
	int rc = sqlite3_open(name.c_str(), &conn2);
	char *zErrMsg = 0;
	sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
	sql = "delete from sequence_profile_map where profile_id =";
	sql += to_string(todelete[i])+";";
	rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
	sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
	
	sqlite3_exec(conn2, "BEGIN TRANSACTION", NULL, NULL, NULL);
	sql = "delete from profile_alignments where id =";
	sql += to_string(todelete[i])+";";
	rc = sqlite3_exec(conn2, sql.c_str(), 0, 0, 0);
	sqlite3_exec(conn2, "COMMIT TRANSACTION", NULL, NULL, NULL);
    }
    get_profile_alignment_nums(updatedprofsnums);
}

void GeneDB::toggle_updated_all_off(){
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string sql = "update alignments set updated = 0;";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
}

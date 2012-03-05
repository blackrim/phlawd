#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>

#include "genedb.h"
#include "libsqlitewrapped.h"
#include "utils.h"

using namespace std;

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
	query.execute("create table alignments(id INTEGER PRIMARY KEY, filename VARCHAR(255));");
	query.execute("create index alignments_filename on alignments(filename);");
	query.execute("create table sequences(id INTEGER PRIMARY KEY, ncbi_id INTEGER, name VARCHAR(255), sequence TEXT);");
	query.execute("create index sequences_ncbi_id on sequences(ncbi_id);");
	query.execute("create index sequences_name on sequences(name);");
	query.execute("create table user_sequences(id INTEGER PRIMARY KEY, ncbi_id INTEGER, name VARCHAR(255), sequence TEXT);");
	query.execute("create index user_sequences_ncbi_id on user_sequences(ncbi_id);");
	query.execute("create index user_sequences_name on user_sequences(name);");
	query.execute("create table sequence_alignment_map(id INTEGER PRIMARY KEY, sequence_id INTEGER, alignment_id INTEGER, sequence TEXT);");
	query.execute("create index sequence_alignment_map_sequence_id sequence_alignment_map(sequence_id);");
	query.execute("create index sequence_alignment_map_alignment_id sequence_alignment_map(alignment_id);");
	query.execute("create table alignment_alignment_map(id INTEGER PRIMARY KEY, parent_id INTEGER, child_id INTEGER);");
	query.execute("create index alignment_alignment_map_parent_id alignment_alignment_map(parent_id);");
	query.execute("create index alignment_alignment_map_child_id alignment_alignment_map(child_id);");
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
	string sql = "insert into user_sequences (ncbi_id,name,sequence) values (";
	sql += user_seqs->at(i).get_comment()+",'";
	sql += user_seqs->at(i).get_id()+"','";
	sql += user_seqs->at(i).get_sequence()+"');";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
	int sid = sqlite3_last_insert_rowid(conn);
	user_seqs->at(i).set_sqlite_id(sid);
    }
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);
}

void GeneDB::add_seqs_to_db(vector<DBSeq> * keep_seqs){
    cout << "adding sequences to database: " << name << endl;
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    for(int i=0;i<keep_seqs->size();i++){
	string sql = "insert into sequences (ncbi_id,name,sequence) values (";
	sql += keep_seqs->at(i).get_ncbi_taxid()+",'";
	sql += keep_seqs->at(i).get_edited_name()+"','";
	sql += keep_seqs->at(i).get_sequence()+"');";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
	int sid = sqlite3_last_insert_rowid(conn);
	keep_seqs->at(i).set_sqlite_id(sid);
    }
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);
}

//TODO: add the reading of the actual alignment file
int GeneDB::add_alignment(string filename, vector<DBSeq> * dbseqs, vector<Sequence> * userseqs){
    cout << "adding alignment to database: " << name << endl;
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string sql = "insert into alignments (filename) values ('";
    sql += filename+"');";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    int alignid = sqlite3_last_insert_rowid(conn);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    std::stringstream ss;
    ss << alignid;
    string alignids = ss.str();
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    for(int i=0;i<dbseqs->size();i++){
	string sql = "insert into sequence_alignment_map (sequence_id, alignment_id,sequence) values (";
	std::stringstream ss2;
	ss2 << dbseqs->at(i).get_sqlite_id();
	string dbsid = ss2.str();
	sql += dbsid+",";
	sql += alignids+",'";
	sql += dbseqs->at(i).get_sequence()+"');";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    }
    for(int i=0;i<userseqs->size();i++){
	string sql = "insert into sequence_alignment_map (sequence_id, alignment_id,sequence) values (";
	std::stringstream ss2;
	ss2 << userseqs->at(i).get_sqlite_id();
	string dbsid = ss2.str();
	sql += dbsid+",";
	sql += alignids+",'";
	sql += userseqs->at(i).get_sequence()+"');";
	rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    }
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);
}

void GeneDB::remove_alignment(int alignid){
    sqlite3 *conn;
    int rc = sqlite3_open(name.c_str(), &conn);
    char *zErrMsg = 0;
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    string sql = "delete from sequence_alignment_map where alignment_id =";
    std::stringstream ss;
    ss << alignid;
    string alignids = ss.str();
    sql += alignids+");";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);

    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    sql = "delete from alignments where id =";
    sql += alignids+");";
    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
}



//void GeneDB::add_profile_alignment(, int child_id1, int child_id2){

//}

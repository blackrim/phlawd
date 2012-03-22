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
 * SQLiteDBController.h
 */

#include <string>
#include <map>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <cstring>
#include <sstream>

using namespace std;


#include "libsqlitewrapped.h"
#include <sqlite3.h>

#include "SQLiteDBController.h"
#include "utils.h"

SQLiteDBController::SQLiteDBController(string dbn):db_name(dbn),division(""),
		count(0){
}

template <class T>
inline std::string to_string (const T& t){
    std::stringstream ss;
    ss << t;
    return ss.str();
}

bool SQLiteDBController::initiate(){
    bool ret = true;
    //check to see if the database exists
    ifstream ifile(db_name.c_str());
    if(ifile){
	cout << "the file: "+db_name+" seems to exist" << endl;
	return false;
    }
    Database conn(db_name);
    Query query(conn);
    query.get_result("create table taxonomy (id INTEGER PRIMARY KEY,ncbi_id INTEGER,name VARCHAR(255),name_class VARCHAR(32),node_rank VARCHAR(32),parent_ncbi_id INTEGER,edited_name VARCHAR(255),left_value INTEGER,right_value INTEGER);" );
    query.free_result();
    query.get_result("CREATE INDEX taxonomy_left_value on taxonomy(left_value);");
    query.free_result();
    query.get_result("CREATE INDEX taxonomy_name on taxonomy(name);");
    query.free_result();
    query.get_result("CREATE INDEX taxonomy_ncbi_id on taxonomy(ncbi_id);");
    query.free_result();
    query.get_result("CREATE INDEX taxonomy_parent_ncbi_id on taxonomy(parent_ncbi_id);");
    query.free_result();
    query.get_result("CREATE INDEX taxonomy_right_value on taxonomy(right_value);");
    query.free_result();

    query.get_result("create table sequence (id INTEGER PRIMARY KEY,ncbi_id INTEGER,accession_id VARCHAR(128),identifier VARCHAR(40),description TEXT,seq LONGTEXT);");
    query.free_result();
    query.get_result("CREATE INDEX sequence_ncbi_id on sequence(ncbi_id);");
    query.free_result();
    query.get_result("CREATE INDEX sequence_accession_id on sequence(accession_id);");
    query.free_result();
    query.get_result("CREATE INDEX sequence_identifier on sequence(identifier);");
    query.free_result();
    return ret;
}

static int callback(void *NotUsed, int argc, char **argv, char **azColName){
    int i;
    for(i=0; i<argc; i++){
	printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
    return 0;
}

void SQLiteDBController::load_seqs(string div,bool downl){
    cout << "loading taxonomy" << endl;
    if (downl == true){
	const char * cmd = "wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";
	cout << "downloading with wget" << endl;
	FILE *fp = popen(cmd, "r" );
	char buff[1000];
	while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	    string line(buff);
	}
	pclose( fp );
	cmd = "tar -xzvf taxdump.tar.gz";
	cout << "untaring" << endl;
	fp = popen(cmd, "r" );
	while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	    string line(buff);
	}
	pclose( fp );
    }
    //read the nodes.dmp
    map<string,string> rank;
    map<string,string> parent_id;
    ifstream infile("nodes.dmp",ios::in);
    string line;
    vector<string> tokens;
    while(getline(infile,line)){
	string del("|");
	tokens.clear();
	Tokenize(line, tokens, del);
	if(tokens.size() > 1){
	    for(int i=0;i<tokens.size();i++){
		TrimSpaces(tokens[i]);
	    }
	    string ncbi_id = tokens[0];
	    rank[ncbi_id] = tokens[2];
	    parent_id[ncbi_id] = tokens[1];
	}
    }
    infile.close();
    //read the names.dmp
    ifstream infile2 ("names.dmp",ios::in);
    sqlite3 *conn;
    int rc = sqlite3_open(db_name.c_str(), &conn);
    char *zErrMsg = 0;

    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    count = 0;
    while(getline(infile2,line)){
	if (count % 100000 == 0){
	    cout << count << endl;
	}
	string del("|");
	tokens.clear();
	Tokenize(line, tokens, del);
	if(tokens.size() > 1){
	    for(int i=0;i<tokens.size();i++){
		TrimSpaces(tokens[i]);
	    }
	    string gin = tokens[0];
	    string nm = tokens[1];//need to double quote the single quotes and backslash the quotes
	    string nm_class = tokens[3];
	    string ednm = tokens[1];//need to edit the names
	    string sql = "insert into taxonomy (ncbi_id,name,name_class,node_rank,parent_ncbi_id,edited_name) values (";
	    sql += gin+",'";
	    sql += nm+"','";
	    sql += nm_class+"','";
	    sql += rank[gin]+"',";
	    sql += parent_id[gin]+",'";
	    sql += ednm+"');";
	    //query.execute(sql);
	    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
	}
	count += 1;
    }
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    infile2.close();
    sqlite3_close(conn);

    cout << "updating left/right values" << endl;
    Database cppconn(db_name);
    Query query(cppconn);
    string cmd = "select ncbi_id,parent_ncbi_id from taxonomy where name_class = 'scientific name';";
    query.get_result(cmd);
    while(query.fetch_row()){
	int nc = query.getval();
	int pc = query.getval();
	if(parent_ncbi_map.count(pc) > 0){
	    parent_ncbi_map[pc].push_back(nc);
	}else{
	    vector<int> tv;
	    tv.push_back(nc);
	    parent_ncbi_map[pc] = tv;

	}
    }
    //get the root and send to rebuild
    count = 0;
    rc = sqlite3_open(db_name.c_str(), &conn);
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    rebuild_tree(1,1,conn);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);

    cout << "loading seqs" << endl;
    division = div;
    if (downl == true){
	string cmd = "wget ftp://ftp.ncbi.nih.gov/genbank/gb"+division+"*.seq.gz";
	cout << "downloading with wget" << endl;
	FILE *fp = popen(cmd.c_str(), "r" );
	char buff[1000];
	while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	    string line(buff);
	}
	pclose( fp );
	cmd = "gunzip -d gb"+division+"*.seq.gz";
	cout << "uncompressing" << endl;
	fp = popen(cmd.c_str(), "r" );
	while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
	    string line(buff);
	}
	pclose( fp );
    }
    vector<string> file_names;
    cout << "getting file names" << endl;
    getdir(".",file_names);
    for(int i=0;i<file_names.size();i++){
	if(file_names[i].find("gb"+div) != string::npos && file_names[i].find(".seq") != string::npos){
	    string filen = file_names[i];
	    cout << filen << endl;
	}
    }
    /*handle = open(filen,"rU")
      con = sqlite3.connect(database)
      curup = con.cursor()
      for i in SeqIO.parse(handle,"gb"):
      acc = i.id
      iden = i.annotations['gi']
      seq = i.seq.tostring()
      desc = i.description
      ncbi_id = ""
      a = i.features[0].qualifiers['db_xref']
      for j in a:
      if 'taxon' in j:
      ncbi_id = j[6:]
      curup.execute("insert into sequence (ncbi_id,accession_id,identifier,description,seq) values (?,?,?,?,?);",(ncbi_id,acc,iden,desc,seq))
      con.commit()*/

    cout << "finished loading" << endl;
}

//private
int SQLiteDBController::rebuild_tree(int gid,int lft,sqlite3 * conn){
    int rgt = lft + 1;
    vector<int>res;
    if (parent_ncbi_map.count(gid) > 0){
	res = parent_ncbi_map[gid];
    }
    for (int i = 0;i < res.size();i++){
	if (res[i] == gid)
	    continue;
	else
	    rgt = rebuild_tree(res[i],rgt,conn);
    }
    string updcmd = "update taxonomy set left_value = "+to_string(lft)+", right_value = "+to_string(rgt)+" where ncbi_id = "+to_string(gid)+";"; //and name_class = scientific name
    //cout << updcmd << endl;
    //exit(0);
    int rc = sqlite3_exec(conn, updcmd.c_str(), 0, 0, 0);
    if (count % 100000 == 0)
	cout << count << endl;
    count += 1;
    return rgt + 1;
}

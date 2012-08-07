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
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <cstring>
#include <sstream>
#include <stdio.h>
#include <regex.h>
#include <unistd.h>

using namespace std;


#include "libsqlitewrapped.h"
#include <sqlite3.h>
#include "GenBankReader.h"
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
    query.get_result("CREATE INDEX taxonomy_edited_name on taxonomy(edited_name);");
    query.free_result();

    query.get_result("create table sequence (id INTEGER PRIMARY KEY,ncbi_id INTEGER,accession_id VARCHAR(128),identifier VARCHAR(40),description TEXT,seq LONGTEXT);");
    query.free_result();
    query.get_result("CREATE INDEX sequence_ncbi_id on sequence(ncbi_id);");
    query.free_result();
    query.get_result("CREATE INDEX sequence_accession_id on sequence(accession_id);");
    query.free_result();
    query.get_result("CREATE INDEX sequence_identifier on sequence(identifier);");
    query.free_result();

    query.get_result("create table information (id INTEGER PRIMARY KEY, name VARCHAR(128), value VARCHAR(128));");
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

string SQLiteDBController::create_name(string & tfilen){
    size_t found;
    found = tfilen.find("\"");
    while(found!=string::npos){
	tfilen.replace(found,1,"'");
	found = tfilen.find("\"",found+1);
    }
    return tfilen;
}

string SQLiteDBController::create_edited_name(string & tfilen){
    size_t found;
    found = tfilen.find(" ");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find(" ",found+2);
    }
    //(take out the parenthetical stuff too)
    found = tfilen.find("(");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find("(",found+1);
    }
    found = tfilen.find(")");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find(")",found+1);
    }
    found = tfilen.find("\"");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find("\"",found+1);
    }
    found = tfilen.find("'");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find("'",found+1);
    }
    found = tfilen.find(".");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find(".",found+1);
    }
    found = tfilen.find("&");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find("&",found+1);
    }
    found = tfilen.find(",");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find(",",found+1);
    }
    found = tfilen.find("\\");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find("\\",found+1);
    }
    found = tfilen.find("/");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find("/",found+1);
    }
    found = tfilen.find(";");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find(";",found+1);
    }
    return tfilen;
}

/*ednm = spls[1].replace('\"',"_").strip()
		ednm = ednm.replace("\'","_")
		ednm = ednm.replace("\\","_")
		ednm = ednm.replace("/","_")
		ednm = ednm.replace("(","_")
		ednm = ednm.replace(")","_")
		ednm = ednm.replace(".","_")
		ednm = ednm.replace("&","_")
		ednm = ednm.replace(",","_")
		ednm = ednm.replace(" ","_")
*/
void SQLiteDBController::load_seqs(string div,bool downl){
    cout << "loading taxonomy" << endl;
    if (downl == true){
	const char * cmd = "wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";
	cout << "downloading with wget" << endl;
	system(cmd);
	cmd = "tar -xzvf taxdump.tar.gz";
	cout << "untaring" << endl;
	system(cmd);
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
	    string nm = create_name(tokens[1]);//need to double quote the single quotes and backslash the quotes
	    string nm_class = tokens[3];
	    string ednm = create_edited_name(tokens[1]);//need to edit the names
	    string sql = "insert into taxonomy (ncbi_id,name,name_class,node_rank,parent_ncbi_id,edited_name) values (";
	    sql += gin+",\"";
	    sql += nm+"\",'";
	    sql += nm_class+"','";
	    sql += rank[gin]+"',";
	    sql += parent_id[gin]+",'";
	    sql += ednm+"');";
	    //query.execute(sql);
	    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
	    //uncomment to get the names that don't commit, mostly bad quotes
//	    if (rc != 0)
//		cout << sql << endl;
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
    //remove files
    remove("citations.dmp");
    remove("division.dmp");
    remove("gc.prt");
    remove("gencode.dmp");
    remove("readme.txt");

    //get the root and send to rebuild
    count = 0;
    rc = sqlite3_open(db_name.c_str(), &conn);
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    rebuild_tree(1,1,conn);
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);

    cout << "loading seqs" << endl;
    division = div;
    vector<string> runall;
    if(division == "met" || division == "all"){
	runall.push_back("pri");runall.push_back("rod");runall.push_back("mam");runall.push_back("vrt");runall.push_back("inv");
	if (division == "all")
	    runall.push_back("pln");runall.push_back("bct");
    }else{
	runall.push_back(division);
    }
    if (downl == true){
	string cmd;
	for(int i = 0;i<runall.size();i++){
	    cmd = "wget ftp://ftp.ncbi.nih.gov/genbank/gb"+runall[i]+"*.seq.gz";
	    cout << "downloading with wget" << endl;
	    system(cmd.c_str());
	    cmd = "gunzip -d gb"+runall[i]+"*.seq.gz";
	    cout << "uncompressing" << endl;
	    system(cmd.c_str());
	}
    }else{
	for(int i = 0;i<runall.size();i++){
	    cmd = "gunzip -d gb"+runall[i]+"*.seq.gz";
	    cout << "uncompressing" << endl;
	    system(cmd.c_str());
	}
    }
    vector<string> file_names;
    cout << "getting file names" << endl;
    getdir(".",file_names);
    for(int i=0;i<file_names.size();i++){
	for(int j=0;j<runall.size();j++){
	    if(file_names[i].find("gb"+runall[j]) != string::npos && file_names[i].substr(file_names[i].size()-4,4)==".seq"){
		string filen = file_names[i];
		cout << filen << endl;
		GenBankReader gbr;
		gbr.parse_file(filen,db_name);
		remove(filen.c_str());
	    }
	}
    }
    
    cout << "merging old names with new names" << endl;
    ifstream infile3 ("merged.dmp",ios::in);
    rc = sqlite3_open(db_name.c_str(), &conn);
    sqlite3_exec(conn, "BEGIN TRANSACTION", NULL, NULL, NULL);
    count = 0;
    while(getline(infile3,line)){
	if (count % 100000 == 0){
	    cout << count << endl;
	}
	string del("|");
	tokens.clear();
	Tokenize(line, tokens, del);
	for(int i = 0;i<tokens.size();i++){TrimSpaces(tokens[i]);}
	if(tokens.size() > 1){
	    string sql = "update sequence set ncbi_id = ";
	    sql += tokens[1];
	    sql += " where ncbi_id = ";
	    sql += tokens[0];
	    sql += ";";
	    rc = sqlite3_exec(conn, sql.c_str(), 0, 0, 0);
	    if (rc != 0)
		cout << sql << endl;
	}
	count += 1;
    }
    infile3.close();
    sqlite3_exec(conn, "COMMIT TRANSACTION", NULL, NULL, NULL);
    sqlite3_close(conn);

    remove("merged.dmp");
    remove("names.dmp");
    remove("nodes.dmp");
    remove("delnodes.dmp");
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

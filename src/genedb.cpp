#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "genedb.h"
#include "libsqlitewrapped.h"

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
    }else{//initializing the database with the SQL structure
	Database conn(name);
	cout << "connected to " << name << endl;
	Query query(conn);
	query.execute("create table alignments(id INTEGER PRIMARY KEY, filename VARCHAR(255));");
	query.execute("create index alignments_filename on alignments(filename);");
	query.execute("create table sequences(id INTEGER PRIMARY KEY, ncbi_id INTEGER, name VARCHAR(255), sequence TEXT);");
	query.execute("create index sequences_ncbi_id on sequences(ncbi_id);");
	query.execute("create index sequences_name on sequences(name);");
	query.execute("create table sequence_alignment_map(id INTEGER PRIMARY KEY, sequence_id INTEGER, alignment_id INTEGER);");
	query.execute("create index sequence_alignment_map_sequence_id sequence_alignment_map(sequence_id);");
	query.execute("create index sequence_alignment_map_alignment_id sequence_alignment_map(alignment_id);");
	query.execute("create table alignment_alignment_map(id INTEGER PRIMARY KEY, parent_id INTEGER, child_id INTEGER);");
	query.execute("create index alignment_alignment_map_parent_id alignment_alignment_map(parent_id);");
	query.execute("create index alignment_alignment_map_child_id alignment_alignment_map(child_id);");
	cout << "database created"<<endl;
    }
}

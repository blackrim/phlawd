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
 * SQLTreeNameConvertor.cpp
 */

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

#include "tree.h"
#include "node.h"
#include "tree_reader.h"

#include "libsqlitewrapped.h"

#include "SQLiteTreeNameConvertor.h"

SQLiteTreeNameConvertor::SQLiteTreeNameConvertor(string filen,string dbs){
	filename = filen;
	db = dbs;
}

Tree * SQLiteTreeNameConvertor::convert(){
	Database conn(db);
	TreeReader nw;
	ifstream infile2(filename.c_str());
	vector<string> lines;
	string line;
	while (getline(infile2, line)){
		lines.push_back(line);
	}
	infile2.close();
	intree = nw.readTree(lines[0]);
	//try {
		cout << "Tree has " << intree->getExternalNodeCount() << " leaves." << endl;
		for(int i=0;i<intree->getExternalNodeCount();i++){
			string dbsearchid = intree->getExternalNode(i)->getName();
			Query query(conn);
			query.get_result("SELECT name FROM taxonomy WHERE ncbi_id = '"+dbsearchid+"' AND name_class = 'scientific name'");
//			StoreQueryResult R = query.store();
			string sname;
			while(query.fetch_row()){
				sname = query.getstr();
			}
			intree->getExternalNode(i)->setName(sname);
		}
	//} catch (Exception e) {
	//	cout << "Error when reading tree." << endl;
	//}
	return intree;
}

void SQLiteTreeNameConvertor::writetree(string outfilename){
	ofstream outfile;
	outfile.open(outfilename.c_str(),ios::out);
	outfile << intree->getRoot()->getNewick(true) << ";" << endl;
	outfile.close();
}

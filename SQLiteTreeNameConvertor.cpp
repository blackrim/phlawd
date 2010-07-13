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
#include <vector>
using namespace std;


#include <Phyl/TreeTemplate.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
using namespace bpp;

#include "libsqlitewrapped.h"

#include "SQLiteTreeNameConvertor.h"

SQLiteTreeNameConvertor::SQLiteTreeNameConvertor(string filen,string dbs){
	filename = filen;
	db = dbs;
}

Tree * SQLiteTreeNameConvertor::convert(){
	Database conn(db);
	Newick * newickReader = new Newick(false); //No comment allowed!
	//try {
		intree = newickReader->read(filename); // Tree in file MyTestTree.dnd
		cout << "Tree has " << intree->getNumberOfLeaves() << " leaves." << endl;
		vector<int> ids = intree->getLeavesId();
		for(int i=0;i<ids.size();i++){
			string dbsearchid = intree->getNodeName(ids[i]);
			Query query(conn);
			query.get_result("SELECT name FROM taxonomy WHERE ncbi_id = '"+dbsearchid+"' AND name_class = 'scientific name'");
//			StoreQueryResult R = query.store();
			string sname;
			while(query.fetch_row()){
				sname = query.getstr();
			}
			intree->setNodeName(ids[i],sname);
		}
	//} catch (Exception e) {
	//	cout << "Error when reading tree." << endl;
	//}
	delete newickReader;
	return intree;
}

void SQLiteTreeNameConvertor::writetree(string outfilename){
	Newick * newickWriter = new Newick(false);
	newickWriter->write(*intree,outfilename,true);
	delete newickWriter;
	delete intree;
}

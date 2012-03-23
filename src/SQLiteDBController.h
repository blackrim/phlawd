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

#ifndef SQLITEDBCONTROLLER_H_
#define SQLITEDBCONTROLLER_H_

#include <string>
#include <vector>
#include <map>

using namespace std;

#include "libsqlitewrapped.h"

class SQLiteDBController{
private:
    string db_name;
    string division;
    int count;
    map<int,vector<int> > parent_ncbi_map;
    int rebuild_tree(int,int,sqlite3 *);
    string create_name(string & tfilen);
    string create_edited_name(string & tfilen);

public:
    SQLiteDBController(string dbn);
    bool initiate();
    void load_seqs(string div,bool downl);
};

#endif /* SQLITEDBCONTROLLER_H_ */

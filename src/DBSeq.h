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
 *  Sequence.h
 */

#ifndef _DBSEQ_H_
#define _DBSEQ_H_

#include <string>

using namespace std;

#include "sequence.h"

class DBSeq:
public Sequence {
private:
    string accession;
    string ncbi_taxid;
    string tax_id;
    string descr;
    string gi;
    string edited_name;

public:
    DBSeq(const string & name, const string & sequence,
	  string acc, string gi, string ncbi_t_id, string t_id, string desc, string ed_name);
    string get_accession();
    string get_gi();
    string get_ncbi_taxid();
    string get_tax_id();
    string get_descr();
    string get_edited_name();
    bool operator==(const DBSeq &other) const;

};

#endif

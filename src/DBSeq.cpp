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
 * DBSeq.cpp
 */

#include <string>

using namespace std;

#include "DBSeq.h"

DBSeq::DBSeq(const string & name, const string & sequence,
	string acc, string gin, string ncbi_t_id, string t_id, string desc)
		: Sequence(name, sequence){
	accession = acc;
	gi = gin;
	ncbi_taxid = ncbi_t_id;
	tax_id = t_id;
	descr = desc;
}

string DBSeq::get_accession(){
	return accession;
}

string DBSeq::get_gi(){
	return gi;
}

string DBSeq::get_ncbi_taxid(){
	return ncbi_taxid;
}

string DBSeq::get_tax_id(){
	return tax_id;
}

string DBSeq::get_descr(){
	return descr;
}

bool DBSeq::operator==(const DBSeq &other) const{
	bool ret = true;
	if(other.accession != accession){
		ret = false;
		return ret;
	}
	if(other.ncbi_taxid != ncbi_taxid){
		ret = false;
		return ret;
	}
	if(other.tax_id != tax_id){
		ret = false;
		return ret;
	}
	if(other.descr != descr){
		ret = false;
		return ret;
	}
	return ret;
}

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

#include <iostream>
#include <cstdio>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#include "SQLiteConstructor.h"
#include "SQLiteProfiler.h"
#include "SQLiteTreeNameConvertor.h"
#include "utils.h"
#include "SQLiteDBController.h"
#include "SmithWatermanGotoh.h"

#include "tree.h"
#include "tree_reader.h"
#include "sequence.h"
#include "omp.h"
#include "SWPS3_matrix.h"

int main(int argc, char* argv[]){
    if(argc != 3){
	cout << "PHLAWD 3.1a" << endl;
	cout << "you need more arguments." << endl;
	cout << "usage: PHLAWD task configfile" << endl;
	cout << "possible tasks include:" << endl;
	cout << "	seqquery     -- just get the distribution of best hits to get cutoffs" << endl;
	cout << "	assemble     -- includes assembling and profiling" << endl;
	cout << "	justprofile  -- just profiles (assumes you have assembled)" << endl;
	cout << "	justassemble -- just assembles (assumes you will profile with justprofile" << endl;
	//cout << "	changenames -- changes from ncbi numbers to names for a file (newick, fasta, newick, phylip)" <<endl;
	cout << "	setupdb      -- tasks performed on the SQLite database" << endl;
	cout << "	outlier      -- analyses on the detection of outliers" << endl;
    }else{
	/*
	 * code to parse through the tasks
	 */
	bool asse = false;
	bool prof = false;
	bool chnames = false;
	bool setupdb = false;
	bool seqquery = false;
	bool outlier = false;
	string argvstr(argv[1]);
	if(argvstr == "assemble"){
	    asse = true;
	    prof = true;
	}else if(argvstr == "justprofile"){
	    prof = true;
	}else if(argvstr == "justassemble"){
	    asse = true;
	}else if(argvstr == "changenames"){
	    chnames = true;
	}else if(argvstr == "setupdb"){
	    setupdb = true;
	}else if(argvstr == "seqquery"){
	    seqquery = true;
	    asse = true;
	}else if(argvstr == "outlier"){
	    outlier = true;
	}else{
	    cout << "you entered an option that doesn't exist: " << argvstr << endl;
	    cout << "try one of these: " <<endl;
	    cout << "	seqquery     -- just get the distribution of best hits to get cutoffs" << endl;
	    cout << "	assemble     -- includes assembling and profiling" << endl;
	    cout << "	justprofile  -- just profiles (assumes you have assembled)" << endl;
	    cout << "	justassemble -- just assembles (assumes you will profile with justprofile" << endl;
	    cout << "	setupdb      -- tasks performed on the SQLite database" << endl;
	    cout << "	outlier      -- analyses on the detection of outliers" << endl;
	}

	if(asse == true || prof == true || outlier == true){
	    /*
	     *
	     */
	    string clade;
	    bool automated = false;
	    vector <string> search;
	    string gene;
	    double mad = 0.01;
	    double coverage = 0.6;
	    double identity = 0.4;
	    string dbtype;
	    string db;
	    int numthreads = 1;
	    string knownfile;
	    bool uselistfile = false;
	    string listfile; //should be genus<space>species
	    bool useexcludelistfile = false;
	    string excludelistfile;
	    bool useexclude_gi_listfile = false;
	    string exclude_gi_listfile;
	    bool useinclude_gi_listfile = false;
	    string include_gi_listfile;
	    bool containshigher = false;
	    bool containswild = false;
	    bool useITS = false;
	    bool updateDB = false;
	    bool updateFILE = false;
	    string updatef = "";
	    string maskurl = "";
	    string genedb;
	    bool usertree = false;//guide tree
	    string usertreefile = "";//guide tree
	    Tree * usertreeobj;
	    bool userfasta = false;
	    string userfastafile = "";
	    bool userskipdb = false;
	    bool userskipsearch = false;
	    //for outlier
	    double taxcutoff = 4;//number above
	    double blcutoff = 100;//proportion above the mean that is considered outlier so 100 times the mean is default
	    bool outliertreerooted = false;
	    string outliertreefile = "";//set below if here
	    bool outliertreecalculated = false;
	    //read file
	    ifstream ifs(argv[2]);
	    string line;
	    while(getline(ifs,line)){
		if (line.size() < 2){
		    continue;
		}
		vector<string> tokens;
		string del("=");
		tokens.clear();
		Tokenize(line, tokens, del);
		for(int j=0;j<tokens.size();j++){
		    TrimSpaces(tokens[j]);
		}
		if(!strcmp(tokens[0].c_str(), "clade")){
		    clade = tokens[1];
		}else if(!strcmp(tokens[0].c_str(), "automated")){
		    automated = true;
		}else if(!strcmp(tokens[0].c_str(), "search")){
		    vector<string> searchtokens;
		    Tokenize(tokens[1], searchtokens, ",");
		    for(int j=0;j<searchtokens.size();j++){
			TrimSpaces(searchtokens[j]);
		    }
		    search = searchtokens;//change to a vector and parse with commas
		}else if(!strcmp(tokens[0].c_str(), "gene")){
		    gene = tokens[1];
		    genedb = tokens[1]+".db";
		}else if(!strcmp(tokens[0].c_str(), "mad")){
		    mad = atof(tokens[1].c_str());
		}else if(!strcmp(tokens[0].c_str(), "coverage")){
		    coverage = atof(tokens[1].c_str());
		}else if(!strcmp(tokens[0].c_str(), "identity")){
		    identity = atof(tokens[1].c_str());
		}else if(!strcmp(tokens[0].c_str(), "db")){
		    db = tokens[1];
		}else if(!strcmp(tokens[0].c_str(), "knownfile")){
		    knownfile = tokens[1];
		}else if(!strcmp(tokens[0].c_str(), "containshigher")){
		    containshigher = true;
		}else if(!strcmp(tokens[0].c_str(), "containswild")){
		    containswild = true;
		}else if(!strcmp(tokens[0].c_str(), "listfile")){
		    listfile = tokens[1];
		    uselistfile = true;
		}else if(!strcmp(tokens[0].c_str(),  "excludelistfile")){
		    excludelistfile = tokens[1];
		    useexcludelistfile = true;
		}else if(!strcmp(tokens[0].c_str(),  "excludegilistfile")){
		    exclude_gi_listfile = tokens[1];
		    useexclude_gi_listfile = true;
		}else if(!strcmp(tokens[0].c_str(),  "includegilistfile")){
		    include_gi_listfile = tokens[1];
		    useinclude_gi_listfile = true;
		}else if(!strcmp(tokens[0].c_str(),  "ITS")){
		    useITS = true;
		}else if(!strcmp(tokens[0].c_str(),  "numthreads")){
		    numthreads = atoi(tokens[1].c_str());
		}else if(!strcmp(tokens[0].c_str(), "gbmask")){
		    maskurl = tokens[1];
		}else if(!strcmp(tokens[0].c_str(), "updateDB")){
		    updateDB = true;
		    cout << "updateDB" << endl;
		}else if(!strcmp(tokens[0].c_str(), "updateFILE")){
		    updateFILE = true;
		    updatef = tokens[1];
		    cout << "updateFILE" << endl;
		    cout << "updated file " << updatef << endl;
		}else if(!strcmp(tokens[0].c_str(), "userguidetree")){
		    usertree = true;
		    usertreefile = tokens[1];
		    cout << "user guide treefile: "<< usertreefile <<endl;
		}else if(!strcmp(tokens[0].c_str(), "userfasta")){
		    userfasta = true;
		    userfastafile = tokens[1];
		    cout << "user fasta file: "<<userfastafile <<endl;
		}else if(!strcmp(tokens[0].c_str(), "userskipdb")){
		    userskipdb = true;
		    cout << "skipping ncbi database check" << endl;
		}else if(!strcmp(tokens[0].c_str(),"userskipsearch")){
		    userskipsearch = true;
		    cout << "skipping ncbi database search all together" << endl;
		}else if(!strcmp(tokens[0].c_str(), "outliertreefile")){
		    outliertreecalculated = true;
		    outliertreefile = tokens[1];
		    cout << "outlier tree file: "<< outliertreefile <<endl;
		}else if(!strcmp(tokens[0].c_str(), "taxcutoff")){
		    taxcutoff = atof(tokens[1].c_str());
		    cout << "setting taxcutoff: " << taxcutoff << endl;
		}else if(!strcmp(tokens[0].c_str(), "blcutoff")){
		    blcutoff = atof(tokens[1].c_str());
		    cout << "setting blcutoff: "<< blcutoff << endl;
		}else if(!strcmp(tokens[0].c_str(),"outliertreerooted")){
		    outliertreerooted = true;
		    cout << "the outlier tree is rooted" << endl;
		}
	    }
	    ifs.close();
	    /*
	     * checks before moving forward
	     */
	    //checking database
	    ifstream ifile(db.c_str());
	    if(!ifile){
		cerr << "genbank database file: " << db << " doesn't exist" << endl;
		exit(0);
	    }
	    //moving on
	    if(asse == true){
		cout << "assembly" << endl;
		SQLiteConstructor * a;
		a = new SQLiteConstructor(clade, search,gene,genedb,mad,coverage,identity,db,knownfile,useITS, numthreads,automated,updateDB,updatef);
		cout << "number of threads: " << a->get_numthreads() << endl;
		omp_set_num_threads(a->get_numthreads());
		cout << "clade name: " << a->get_cladename() << endl;
		for(int i=0;i<a->get_search().size(); i++){
		    cout << "search: " << a->get_search()[i] << endl;
		}
		cout << "gene name: " << a->get_genename() << endl;
		cout << "gene database: " << a->get_genedb() << endl;
		cout << "mad cutoff: " << a->get_madcutoff() << endl;
		cout << "coverage: " << a->get_coverage() << endl;
		cout << "identity: " << a->get_identity() << endl;
		if(uselistfile == true){
		    cout << "using names in list file: " << listfile << endl;
		    a->set_only_names_from_file(listfile,containshigher,containswild);
		}
		if(useexcludelistfile == true){
		    cout << "excluding names in list file: " << excludelistfile << endl;
		    a->set_exclude_names_from_file(excludelistfile);
		}
		if(useexclude_gi_listfile == true){
		    cout << "excluding gi's in list file: " << exclude_gi_listfile << endl;
		    a->set_exclude_gi_from_file(exclude_gi_listfile);
		}
		if(useinclude_gi_listfile == true){
		    cout << "including gi's in list file: " << include_gi_listfile << endl;
		    a->set_include_gi_from_file(include_gi_listfile);
		}
		if(maskurl.size() > 0){
		    cout << "getting genbank mask from: " << maskurl << endl;
		    //vector<string> query_gis = query_mask(maskurl);
		}
		if(useITS == true){
		    cout << "using ITS mode: true" << endl;
		    cout << "warning: highly experimental" << endl;
		}
		if(usertree == true){
		    cout << "using user guide tree: "<< usertreefile <<endl;
		    a->set_user_guide_tree(usertreefile,userskipdb);
		}
		if(userfasta == true){
		    cout << "using user fasta file: "<< userfastafile << endl;
		    a->set_user_fasta_file(userfastafile,userskipdb);
		}
		if(userskipsearch == true){
		    a->set_user_skip_search();
		}
		if(seqquery == true){
		    a->set_justseqquery(true);
		}
		int ret = a->run();
		if (ret == 2){
		    updateDB = false; 
		    a->run();
		}
		if(usertree == true){
		    usertreeobj = a->get_user_guide_tree_obj();
		    if (usertreeobj == NULL)//can happen if little overlap between tree and dataset
			usertree = false;
		}
		updateDB = a->get_updatestatus(); //if there was an enormous update, it will make sure to note that
		delete(a);
	    }
	    /*
	     * Profiling
	     */
	    if(prof == true){
		SQLiteProfiler * b;
		b = new SQLiteProfiler(gene,genedb,clade,db,automated,updateDB);
		b->prelimalign();
		if(usertree == true && asse == true){
		    b->set_user_guide_tree(usertreeobj);
		}
		b->run();
		delete(b);
	    }
	    /*
	     * change the file (tree or data) from id to names
	     */

	    /*
	     * this needs a lot of editing -- need multiple trees, fasta files
	     */
	    if(outlier == true){
		cout << "attempting outlier detection" << endl;
		//use the major clade group to pick an outgroup which would just be a random child that has taxa
		string treefile;
		if(outliertreecalculated == false){
		    cout << "converting "<< gene <<".FINAL.aln to phylip format for tree building"<<   endl;
		    convert_to_phylip((gene+".FINAL.aln"),(gene+".FINAL.phy"));
		    cout << "constructing a tree for detection" << endl;
		    remove(("RAxML_info."+gene+".outlierdet").c_str());
		    string treerunstr = "raxmlHPC-PTHREADS-SSE3 -T "+int_to_string(numthreads)+" -p 12345 -f E -m GTRCAT -c 8 -s "+(gene+".FINAL.phy")+" -n "+gene+".outlierdet";
		    system(treerunstr.c_str());
		    treefile = "RAxML_bestTree."+gene+".outlierdet";
		}else{
		    treefile = outliertreefile;
		}
		TreeReader nw;
		cout << "outlier reading tree file" << endl;
		ifstream tinfile(treefile.c_str());
		if (!tinfile){
		    cerr <<"outlier tree: "<<treefile<< " doesn't exist"<<endl;
		    exit(0);
		}
		string line;
		getline(tinfile,line);
		tinfile.close();
		Tree * outliertree = nw.readTree(line);
		if (outliertreerooted == false){
		    cout << "attempting to root based on taxonomy" << endl;	
		    get_earliest_branch_representation(db,clade,outliertree);
		}
		cout << "calculating taxonomic outliers" << endl;
		get_taxonomic_outliers(outliertree,db,taxcutoff,gene);
		cout << "calculating branch length outliers" << endl;
		get_branch_length_outliers(outliertree,blcutoff,gene);
		cout << "generating outfile with outliers removed" << endl;
	    }
	}else if (chnames == true){//change names == true
	    cout << "changing names in tree" << endl;
	    string dbtype;
	    string db;
	    string infile;
	    string outfile;
	    //read file
	    ifstream ifs(argv[2]);
	    string line;
	    while(getline(ifs,line)){
		vector<string> tokens;
		string del("=");
		tokens.clear();
		Tokenize(line, tokens, del);
		for(int j=0;j<tokens.size();j++){
		    TrimSpaces(tokens[j]);
		}
		if(!strcmp(tokens[0].c_str(), "infile")){
		    infile = tokens[1];
		}else if(!strcmp(tokens[0].c_str(), "outfile")){
		    outfile = tokens[1];
		}else if(!strcmp(tokens[0].c_str(),  "db")){
		    db = tokens[1].c_str();
		}
	    }
	    ifs.close();

	    //mysql == true
	    cout << "using sqlite" << endl;
	    SQLiteTreeNameConvertor * c;
	    c = new SQLiteTreeNameConvertor(infile,db);
	    c->convert();
	    c->writetree(outfile);
	    delete c;
	}else if(setupdb == true){
	    cout << "setting up database" << endl;
	    string dbname;
	    string division = "";
	    bool download = false;
	    //read file
	    ifstream ifs(argv[2]);
	    string line;
	    while(getline(ifs,line)){
		vector<string> tokens;
		string del("=");
		tokens.clear();
		Tokenize(line, tokens, del);
		for(int j=0;j<tokens.size();j++){
		    TrimSpaces(tokens[j]);
		}
		if(!strcmp(tokens[0].c_str(), "db")){
		    dbname = tokens[1];
		}else if(!strcmp(tokens[0].c_str(), "division")){
		    division = tokens[1];
		}else if(!strcmp(tokens[0].c_str(),  "download")){
		    download = true;
		}
	    }
	    ifs.close();

	    SQLiteDBController * c;
	    c = new SQLiteDBController(dbname);
	    bool ret = c->initiate();
	    if (ret == false){
		cout << "database exists: exiting" << endl;
		exit(0);
	    }
	    if(division.size () > 1){
		c->load_seqs(division,download);
	    }else{
		cout << "you need to add a division=<division,like pln> to "<< argv[2] <<endl;
		exit(0);
	    }
	    //need to add updating

	    delete c;
	}
    }
}



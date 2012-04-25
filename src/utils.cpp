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
 *  utils.cpp
 */

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include "sequence.h" 
#include "fasta_util.h"
#include <fstream>
#include "libsqlitewrapped.h"
#include "tree.h"

using namespace std;

#include "SWPS3_matrix.h"

template <class T>
inline std::string to_string (const T& t){
    std::stringstream ss;
    ss << t;
    return ss.str();
}


void Tokenize(const string& str, vector<string>& tokens,
                      const string& delimiters = " "){
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void TrimSpaces( string& str)  {
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if(( string::npos == startpos ) || ( string::npos == endpos))
    {
	str = "";
    }
    else
	str = str.substr( startpos, endpos-startpos+1 );
    /*
    // Code for  Trim Leading Spaces only
    size_t startpos = str.find_first_not_of(� \t�); // Find the first character position after excluding leading blank spaces
    if( string::npos != startpos )
    str = str.substr( startpos );
    */
    /*
    // Code for Trim trailing Spaces only
    size_t endpos = str.find_last_not_of(� \t�); // Find the first character position from reverse af
    if( string::npos != endpos )
    str = str.substr( 0, endpos+1 );
    */
}

double median(vector<double> x){
    sort(x.begin(),x.end());
    double n = x.size();
    return ( x [ n / 2.0 ] + x [ n / 2.0 - 1.0] ) / 2.0F;
}

double mean(vector<double> & x){
    double size = x.size();
    double ret = 0;
    for(int i=0;i<x.size();i++){
	ret += x[i]/size;
    }
    return ret;
}

double stdev(vector<double> & x){
    double mn = mean(x);
    vector <double> devs(x.size());
    for(int i=0;i<x.size();i++){
	devs[i] = x[i] - mn;
	devs[i] = (devs[i]*devs[i]);
    }
    double sum = 0;
    for(int i=0;i<x.size();i++){
	sum += devs[i];
    }
    sum = (sum/(x.size()-1));
    return sqrt(sum);
}

int getdir (string dir, vector<string> &files){
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
    	/*
    	 * no dot files
    	 */
    	bool test = false;
    	if(*(string(dirp->d_name).begin())=='\.')
	    test = true;
    	if(test == false){
	    files.push_back(string(dirp->d_name));
    	}
    }
    closedir(dp);
    return 0;
}

void splitstring(string str, string seperater, string &first, string &second){
    int i = (int)str.find(seperater); //find seperator
    if(i != -1){
	int y = 0;
	if(!str.empty()){
	    while(y != i){
		first += str[y++]; //creating first string
	    }
	    y = y+(int)seperater.length(); //jumping forward seperater length
	    while(y != str.length()){
		second += str[y++]; //creating second string
	    }

	}
    }
    else{
	first = str;
	second = "NULL"; //if seperator is not there then second string == null
    }
}

void fix_bad_chars(string & tfilen){
    size_t found;
    found = tfilen.find(" ");
    while(found!=string::npos){
	tfilen.replace(found,1,"\\ ");
	found = tfilen.find(" ",found+2);
    }
    //(take out the parenthetical stuff too)
    found = tfilen.find("(");
    while(found!=string::npos){
	tfilen.replace(found,1,"\\(");
	found = tfilen.find("(",found+2);
    }
    found = tfilen.find(")");
    while(found!=string::npos){
	tfilen.replace(found,1,"\\)");
	found = tfilen.find(")",found+2);
    }
}

void fix_bad_chars_for_seq_names(string & tfilen){
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
	found = tfilen.find("(",found+2);
    }
    found = tfilen.find(")");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find(")",found+2);
    }
    found = tfilen.find(".");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find(".",found+2);
    }
    found = tfilen.find("&");
    while(found!=string::npos){
	tfilen.replace(found,1,"_");
	found = tfilen.find("&",found+2);
    }
}

string int_to_string(int n){
   stringstream ss;
   ss << n;
   return ss.str();
}


//get the swps3 sequence score
int get_swps3_score_and_rc_cstyle(SBMatrix mat,  Sequence * inseq1, Sequence * inseq2){
    return swps3_maxscores(mat, inseq1,inseq2);
}

/*
 * This will query the mask given the mask url provided in the config file.
 * It relies on curl to be installed and in the path. Though this is a 
 * relatively low bar considering it is installed generally by default in 
 * both MAC and Linux (basically all UNIX) machines.
 *
 * Would rather use libcurl but might be more difficult for source
 * installers.
 *
 * url = the input url
 * output = vector of the gi strings
 */
vector<string> query_mask(string url){
    string systemcall = "curl "+url;
    systemcall += " > gbmask.downloaded";
    cout << systemcall << endl;
    //need a check here to make sure that it is valid
    system(systemcall.c_str());
    vector<string> returngis;

    return returngis;
}

void convert_to_phylip(string infile,string outfile){
    FastaUtil fu;
    vector<Sequence> tseqs;
    fu.readFile(infile,tseqs);
    ofstream ofile(outfile.c_str());
    ofile << tseqs.size()<<"\t" <<tseqs[0].get_sequence().size() << endl;
    for (int i=0;i<tseqs.size();i++){
	ofile << tseqs[i].get_id() << "\t"<<tseqs[i].get_sequence() << endl;
    }
    ofile.close();
    
}

/*
 * for calculating outliers
 * this will reroot based on one of the children from the root id, the smallest
 * that is at least one of many children, requires that the tips of the tree be
 * ncbi_ids and the root be the clade name
 */
void get_earliest_branch_representation(string ncbidb,string root, Tree * tree){
    Database conn(ncbidb);
    Query query(conn);
    string rootid;
    query.get_result(("SELECT ncbi_id FROM taxonomy WHERE name = '"+root+"'").c_str());
    while(query.fetch_row()){
	rootid = to_string(query.getval());
	break;
    }
    cout << root << endl;
    query.free_result();
    vector<string> trncbis;
    vector<int> trleft;
    vector<int> trright;
    for (int i=0;i<tree->getExternalNodeCount();i++){
	trncbis.push_back(tree->getExternalNode(i)->getName());
	string searchstr = "select left_value,right_value from taxonomy where ncbi_id = "+tree->getExternalNode(i)->getName()+" and name_class='scientific name';";
	query.get_result(searchstr);
	while(query.fetch_row()){
	    trleft.push_back(atoi(to_string(query.getval()).c_str()));
	    trright.push_back(atoi(to_string(query.getval()).c_str()));
	}
	query.free_result();
    }
    bool going = true;
    //keep going until have children that have more than one representative
    string curid = rootid;
    vector<int> ncbis;
    vector<int> ncbis_l;
    vector<int> ncbis_r;
    while(going){
	string searchstr = "select left_value,right_value,ncbi_id from taxonomy where parent_ncbi_id = "+curid+" and name_class='scientific name' and ncbi_id > 1;";
	query.get_result(searchstr);
	while(query.fetch_row()){
	    int left_v = atoi(to_string(query.getval()).c_str());
	    int right_v = atoi(to_string(query.getval()).c_str());
	    int ncbiid = atoi(to_string(query.getval()).c_str());
	    bool matched = false;
	    for(int i=0;i<trncbis.size(); i++){
		if(trleft[i] > left_v && trright[i] < right_v){
		    matched = true;
		    break;
		}
	    }
	    if (matched == true){
		ncbis.push_back(ncbiid);
		ncbis_l.push_back(left_v);
		ncbis_r.push_back(right_v);
	    }
	}
	query.free_result();
	if(ncbis.size() > 1){
	    going = false;
	    break;
	}else{//has to be one
	    curid = ncbis[0];
	    ncbis.clear();
	    ncbis_l.clear();
	    ncbis_r.clear();
	}
    }
    //get smallest max distance or first to have distance or size 0/1
    int lowestscore=trncbis.size();
    int lowest;
    vector<Node *> final_node_vec;
    for (int i=0;i<ncbis.size();i++){
	vector<Node *> node_vec;
	for(int j=0;j<tree->getExternalNodeCount();j++){
	    if(trleft[j] > ncbis_l[i] && trright[j] < ncbis_r[i]){
		node_vec.push_back(tree->getExternalNode(j));
	    }
	}
	if (node_vec.size() < lowestscore){
	    if (tree->getMRCA(node_vec)==tree->getRoot()){
		continue;
	    }
	    lowest = i;
	    lowestscore = node_vec.size();
	    final_node_vec = node_vec;
	    if (node_vec.size() == 1)
		break;
	}
    }
    //get the mrca of the final_node_vec and reroot
    cout << "rerooting with " << ncbis[lowest] << endl;
    Node * mrca = tree->getMRCA(final_node_vec);
    tree->reRoot(mrca);
}

//for distances for taxonomic outliers
vector<int> get_left_right_exclude(vector<int> * lefts, vector<int> * rights, vector<int> * exlefts, vector<int> * exrights){
    vector<int> lf_rt;
    lf_rt.push_back(1e+20);
    lf_rt.push_back(0);
    for (int i=0;i<lefts->size();i++){
        if (count(exlefts->begin(), exlefts->end(),lefts->at(i))==0)
	    if (lefts->at(i) < lf_rt[0])
		   lf_rt[0] = lefts->at(i);
    }
    for (int i=0;i<rights->size();i++){
        if (count(exrights->begin(),exrights->end(),rights->at(i))==0)
            if (rights->at(i) > lf_rt[1])
	lf_rt[1] = i;
}
return lf_rt;
}
//for distance for taxonomic outliers
int get_distance_from_child_to_parent(Query * qu, string child, string parent){
    int distance = 0;
    if (child == parent)
        return 0;
    string curid = child;
    while (curid != "1"){
        string sqlsearch = "select parent_ncbi_id from taxonomy where ncbi_id = "+curid+";";
        qu->get_result(sqlsearch);
        string nid;
	while(qu->fetch_row()){
	    nid = to_string(qu->getval());
	}
	qu->free_result();
        distance += 1;
        if (nid == parent)
            break;
        curid = nid;
    }
    return distance;
}

//postorder calculator
void get_distances(Node * curnode,map<Node *,int> * distances,map<Node*,vector<int> > * trleft_right, Query * qu){
    for(int i=0;i<curnode->getChildCount();i++){
	get_distances(curnode->getChild(i),distances,trleft_right,qu);
    }
    if (curnode->getParent() == NULL)
	return;
    vector<string> leaves_to_ex;
    string ilabel;
    if (curnode->getChildCount()==0){
	ilabel = curnode->getName();
	leaves_to_ex.push_back(ilabel);
    }else{
	ilabel = curnode->getName();
	vector<Node *> lvs = curnode->get_leaves();
	for (int j=0;j<lvs.size();j++){
	    string jlabel = lvs[j]->getName();
	    leaves_to_ex.push_back(jlabel);
	}
    }
    vector<int>lefts;
    vector<int>rights;
    vector<int>exlefts;
    vector<int>exrights;
    Node * testnode = curnode->getParent();
    if (curnode->getParent()->get_num_leaves() == 2 || (curnode->getParent()->get_num_leaves()-leaves_to_ex.size() == 1)){
	Node *testnode = curnode->getParent()->getParent();
	if (testnode == NULL)
	    return;
    }
    vector<Node *> tleaves = testnode->get_leaves();
    for(int j=0;j<tleaves.size();j++){
	string testlabel = tleaves[j]->getName();
	if (count(leaves_to_ex.begin(),leaves_to_ex.end(),testlabel)>0){
	    exlefts.push_back((*trleft_right)[tleaves[j]][0]);
	    exrights.push_back((*trleft_right)[tleaves[j]][1]);
	}else{
	    lefts.push_back((*trleft_right)[tleaves[j]][0]);
	    rights.push_back((*trleft_right)[tleaves[j]][1]);
	}
    }
    vector<int> lf_rt = get_left_right_exclude(&lefts,&rights,&exlefts,&exrights);
    string sqlstring = "select ncbi_id from taxonomy where left_value < "+to_string(lf_rt[0])+" and right_value > "+to_string(lf_rt[1])+" LIMIT 1;";
    qu->get_result(sqlstring);
    string commonparent;
    while(qu->fetch_row()){
	commonparent = to_string(qu->getval());
	break;
    }
    qu->free_result();
    int distance = get_distance_from_child_to_parent(qu,commonparent,testnode->getName());
    (*distances)[curnode] = distance;
}

//for taxonomic outliers
void get_suggested_clips(map<Node *, int> * distances, Tree * tree, double meand,map<Node *,int> *marked, map<string, int>* finallvs, double cutoff){
    //dont' know if the preorder is necessary
    //    for i in tree.iternodes(order=0):
    for(int i=0;i<tree->getNodeCount();i++){
        if (tree->getNode(i)->getParent() == NULL || distances->count(tree->getNode(i))==0){
            //cout << i << "not in distances" << endl;
	    continue;
	}
        /* check if the other child is larger
	 * this typically means that one child is driving down the other a little
	 * and the big one is the bad one
	 */
	int testdistance = (*distances)[tree->getNode(i)];
	if ( testdistance > cutoff+meand){
	    Node * othernode = tree->getNode(i)->getSister();
	    if ((*distances)[othernode] > testdistance){
		continue;
	    }else{
		(*marked)[tree->getNode(i)] = testdistance;
		// good to take out parents because of trickle
		Node * curnode = tree->getNode(i);
		while (curnode->getParent() != NULL){
		    curnode = curnode->getParent();
		    if (marked->count(curnode)>0){
			// if marked_nodes[curnode] <= distances[i]:
			map<Node*,int>::iterator it;
			it = marked->find(curnode);
			marked->erase(it);
		    }		
		}		
	    }
	}
    }

    map<Node*,int>::iterator it;
    for(it=marked->begin();it != marked->end();it++){
	vector<Node *> lvs = (*it).first->get_leaves();
	for(int j = 0;j<lvs.size();j++){
	    if (finallvs->count(lvs[j]->getName())==0)
		(*finallvs)[lvs[j]->getName()] = (*marked)[(*it).first];
	    else
		if ((*finallvs)[lvs[j]->getName()] < (*marked)[(*it).first])
		    (*finallvs)[lvs[j]->getName()] = (*marked)[(*it).first];
	}
    }
}

void get_taxonomic_outliers(Tree * tree, string ncbidb, double cutoff,string gene){
    Database conn(ncbidb);
    Query query(conn);
    cout << "getting tree tip left right values" << endl;
    map<Node*,vector<int> > trleft_right;
    for (int i=0;i<tree->getExternalNodeCount();i++){
	string searchstr = "select left_value,right_value from taxonomy where ncbi_id = "+tree->getExternalNode(i)->getName()+" and name_class='scientific name';";
	query.get_result(searchstr);
	while(query.fetch_row()){
	    vector<int> ttemp;
	    ttemp.push_back(atoi(to_string(query.getval()).c_str()));
	    ttemp.push_back(atoi(to_string(query.getval()).c_str()));
	    trleft_right[tree->getExternalNode(i)] = ttemp;
	}
	query.free_result();
    }
    cout << trleft_right.size() << " leaves in the tree" << endl;
    map<string,string> labelnames;
    for(int i=0;i<tree->getInternalNodeCount();i++){
        vector<string> leavesnms;
        int lowestleft = 1e+100;
        int highestright = 0;
	vector<Node *> lf_nodes = tree->getInternalNode(i)->get_leaves();
	for(int j=0;j<lf_nodes.size();j++){
	    leavesnms.push_back(lf_nodes[j]->getName());
	    if(trleft_right[lf_nodes[j]][0] < lowestleft)
		lowestleft = trleft_right[lf_nodes[j]][0];
	    if(trleft_right[lf_nodes[j]][1] > highestright)
		highestright = trleft_right[lf_nodes[j]][1];
	}
        string searchstr = "select ncbi_id,edited_name from taxonomy where left_value < "+to_string(lowestleft)+" and right_value > "+to_string(highestright)+" and name_class = 'scientific name' LIMIT 1;";
	query.get_result(searchstr);
	int commonparent;
	string ed_name;
	while(query.fetch_row()){
	    commonparent = atoi(to_string(query.getval()).c_str());
	    ed_name = to_string(query.getstr());
	    break;
	}
        labelnames[to_string(commonparent)] = ed_name;
	tree->getInternalNode(i)->setName(to_string(commonparent));
    }
    map<Node *,int> distances;
    get_distances(tree->getRoot(),&distances,&trleft_right,&query);
    map<Node *,int>::iterator it;
    double meandistance = 0;
    for (it = distances.begin();it != distances.end(); it ++){
	meandistance += ((*it).second/float(distances.size()));
    }
    cout <<"mean distance: " <<  meandistance << endl;
    cout << "cutoff: " << cutoff+meandistance << endl;
    map<Node*, int> marked;
    map<string,int> finallvs;
    get_suggested_clips(&distances, tree, meandistance,&marked,&finallvs, cutoff);
    string outfilename = gene+".taxoutliers";
    cout << "writing suggested taxonomic outliers to file: "<< outfilename  << endl;
    ofstream toutf;
    toutf.open((outfilename).c_str(),ios::out);
    map<string, int>::iterator it2;
    for(it2=finallvs.begin();it2!=finallvs.end();it2++){
	toutf << (*it2).first << "\t" << (*it2).second << endl;
    }	    
    toutf.close();
    outfilename = outfilename+".tre";
    ofstream toutft;
    toutft.open((outfilename).c_str(),ios::out);
    cout << "plotting tax outliers to tree: " << outfilename << endl;
    toutft << "#NEXUS\nbegin trees;" << endl;
    for (int i=0;i<tree->getNodeCount();i++){
        if (tree->getNode(i)->getParent() == NULL || distances.count(tree->getNode(i)) == 0)
            continue;
	if(distances[tree->getNode(i)]> meandistance+cutoff)
	    tree->getNode(i)->setColor("#-2610389");
	else
	    tree->getNode(i)->setColor("#-13903798");
	if(tree->getNode(i)->getChildCount() > 0){
	    tree->getNode(i)->setName(to_string(distances[tree->getNode(i)]));
	}
    }
    toutft << "\ttree tree1 = [&R] " << tree->getRoot()->getNewickColor() << ";" << endl;
    toutft << "end;" << endl;
    toutft.close();
    for (int i=0;i<tree->getNodeCount();i++){
	tree->getNode(i)->setColor("");
	if(tree->getNode(i)->getChildCount() > 0){
	    tree->getNode(i)->setName("");
	}
    }
}


void get_branch_length_outliers(Tree * tree,double cutoff,string gene){
    double meanbl = 0;
    for (int i=0;i<tree->getNodeCount();i++){
	if (tree->getNode(i)->getParent()!= NULL)
	    meanbl += (tree->getNode(i)->getBL()/float(tree->getNodeCount()));
    }
    int totaltips = tree->getExternalNodeCount();
    cout << "mean bl: " << meanbl << endl;
    cutoff = cutoff*meanbl;
    cout << "cutoff: " << cutoff << endl;
    vector<Node *> clips;
    for (int i=0;i<tree->getNodeCount();i++){
	if (tree->getNode(i)->getParent()!= NULL){
	    if(tree->getNode(i)->getParent()->getParent()!=NULL){
		if(tree->getNode(i)->getBL() > meanbl + cutoff){
		    vector<Node *> lvs = tree->getNode(i)->get_leaves();
		    if (lvs.size()/float(totaltips) > 0.1){
			Node * sis = tree->getNode(i)->getSister();
			if (sis->get_leaves().size() < lvs.size())
			    lvs = sis->get_leaves();
		    }
		    for(int j=0;j<lvs.size();j++){
			if(count(clips.begin(),clips.end(),lvs[j])==0)
			    clips.push_back(lvs[j]);
		    }	
		}
	    }
	}
    }
    string outfilename = gene+".bloutliers";
    cout << "writing bl outliers to file: "<< outfilename << endl;
    ofstream toutft;
    toutft.open((outfilename).c_str(),ios::out);
    vector<string> names;
    for (int i=0;i<clips.size();i++){
	toutft << clips[i]->getName() << endl;
	clips[i]->setColor("#-2610389");
    }
    toutft.close();
    outfilename += ".tre";
    cout << "plotting bl outliers to tree: " << outfilename << endl;
    ofstream toutf;
    toutf.open((outfilename).c_str(),ios::out);
    toutf << "#NEXUS\nbegin trees;" << endl;
    toutf << "\ttree tree1 = [&R] " << tree->getRoot()->getNewickColor() << ";" << endl;
    toutf << "end;" << endl;
    toutf.close();
}

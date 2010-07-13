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
 *  Same_seq_pthread.cpp
 */

#include <Seq/Sequence.h>
#include <Seq/containers>
#include <Seq/ioseq>
#include <Seq/alphabets>
using namespace bpp;

#include "DBSeq.h"
#include "utils.h"

#include <cstdio>
#include <string>
#include <vector>
using namespace std;

#include "Same_seq_pthread.h"

vector<double> get_blast_score_and_rc_cstyle(Sequence inseq1, DBSeq inseq2, bool * prevcomp, int thread){
	const NucleicAlphabet * alphabet = new DNA();
	vector<double> retvalues;
	Fasta seqwriter1;
	Fasta seqwriter2;
	VectorSequenceContainer * sc1 = new VectorSequenceContainer(alphabet);
	sc1->addSequence(inseq1);
	VectorSequenceContainer * sc2 = new VectorSequenceContainer(alphabet);
	sc2->addSequence(inseq2);

	std::string ts;
	std::stringstream out;
	out << thread;
	ts = out.str();
	const string fn1 = "seq1"+ts;
	const string fn2 = "seq2"+ts;
	seqwriter1.write(fn1,*sc1);
	seqwriter2.write(fn2,*sc2);
	delete sc1;
	delete sc2;
	delete alphabet;
//	string cmd = "bl2seq -i seq1 -j seq2 -p blastn -D 1";
	double coverage = 0;
	double identity = 0;
	string line;

	std::string stcmd = "bl2seq -i seq1";
	stcmd += ts;
	stcmd += " -j seq2";
	stcmd += ts;
	stcmd += " -p blastn -D 1";

	//cout <<"thread: "<<thread << " " << stcmd << endl;

	const char * cmd = stcmd.c_str();
	FILE *fp = popen(cmd, "r" );
	char buff[1000];
	vector<string> tokens;
	if(fp){
		while ( fgets( buff, sizeof buff, fp ) != NULL ) {//doesn't exit out
			string line(buff);
			size_t found=line.find("#");
			if (found==string::npos){
				string del("\t");
				Tokenize(line, tokens, del);
				coverage = coverage + (strtod(tokens[3].c_str(),NULL)-strtod(tokens[4].c_str(),NULL));
				double tid = (strtod(tokens[2].c_str(),NULL)/100.0)*
					((strtod(tokens[3].c_str(),NULL)-strtod(tokens[4].c_str(),NULL))/(int)inseq1.toString().length());

				if (tid > identity){
					identity = tid;
					//test
					if (strtod(tokens[8].c_str(),NULL)>strtod(tokens[9].c_str(),NULL)){
						*prevcomp=true;
					}else{
						*prevcomp=false;
					}
				}
				//cout << tid <<" " << identity << " " << *prevcomp << " " << strtod(tokens[8].c_str(),NULL) << " " << strtod(tokens[9].c_str(),NULL) << endl;
				//cout << "try " << line;
				//cout << "tokens " << tokens[2] << " " << tokens[3] << endl;

			}
		}
		pclose( fp );
	}
	//delete cmd;
	//delete fp;
	if (tokens.size() < 1){
		return retvalues;
	}/*else{
		if (strtod(tokens[8].c_str(),NULL)>strtod(tokens[9].c_str(),NULL)){
			*prevcomp=true;
		}else{
			*prevcomp=false;
		}
	}*/
	retvalues.push_back(identity);
	retvalues.push_back(coverage/(int)inseq1.toString().length());
	return retvalues;
	//return (float(maxident/100.0),float(coverage/len(seq1.seq.tostring())),rc)
}

void * Same_seq_pthread_go(void *threadarg){
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	int thread_id;
	vector<DBSeq> seqs;
	vector<DBSeq> keep_seqs;
	vector<bool> keep_rc;
	int reports;
	double coverage;
	double identity;
	OrderedSequenceContainer * known_seqs;

	thread_id = my_data->thread_id;
	seqs = my_data->seqs;
	keep_seqs = my_data->keep_seqs;
	keep_rc = my_data->keep_rc;
	reports = my_data->reports;
	coverage = my_data->coverage;
	identity = my_data->identity;
	known_seqs = my_data->known_seqs;

	for (int i=0;i<seqs.size();i++){
		if(i%reports == 0){
			cout << i << endl;
		}
		double maxide = 0;
		double maxcov = 0;
		bool rc = false;
		for (int j=0;j<known_seqs->getNumberOfSequences();j++){
			bool trc = false;
			//TODO : there was a pointer problem here
			vector<double> ret = get_blast_score_and_rc_cstyle(known_seqs->getSequence(j), seqs[i],&trc, thread_id); //should be pointer?
			if (ret.size() > 1){
				/*if (ret[0] >maxide){
					maxide = ret[0];
				}
				if (ret[1] > maxcov){ // should these be in the same conditional statement
					maxcov = ret[1];
					rc = trc;//need to get it somewhere else -- pointer probably
				}*/
				//cout << ret[0] << " " << ret[1] << endl;
				if (ret[0] >maxide && ret[1] > maxcov){
					maxide = ret[0];
					maxcov = ret[1];
					rc = trc;
				}
			}
		}
		if (maxide >= identity && maxcov >= coverage){
			keep_seqs.push_back(seqs[i]);
			keep_rc.push_back(rc);
			//cout << keep_seqs.size() << endl;
		}
	}
	my_data->keep_seqs = keep_seqs;
	my_data->keep_rc = keep_rc;
}


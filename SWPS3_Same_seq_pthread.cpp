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
#include <Seq/sequences>
using namespace bpp;

#include "DBSeq.h"
#include "utils.h"

#include <cstdio>
#include <string>
#include <vector>
using namespace std;

#include "SWPS3_matrix.h"
#include "SWPS3_Same_seq_pthread.h"


int get_swps3_score_and_rc_cstyle(SBMatrix mat, const Sequence * inseq1, const Sequence * inseq2){

	return swps3_maxscores(mat, inseq1,inseq2);
}

void * SWPS3_Same_seq_pthread_go(void *threadarg){
	struct SWPS3_thread_data *my_data;
	my_data = (struct SWPS3_thread_data *) threadarg;
	int thread_id;
	vector<DBSeq> seqs;
	vector<DBSeq> keep_seqs;
	vector<bool> keep_rc;
	int reports;
	double identity;
	OrderedSequenceContainer * known_seqs;
	thread_id = my_data->thread_id;
	seqs = my_data->seqs;
	keep_seqs = my_data->keep_seqs;
	keep_rc = my_data->keep_rc;
	reports = my_data->reports;
	identity = my_data->identity;
	known_seqs = my_data->known_seqs;

	/*
	 * get the best score for each sequence -- the lowest of known to known
	 */
	vector<int> scores;
	SBMatrix mat = swps3_readSBMatrix( "EDNAFULL" );
	for(int i=0;i<known_seqs->getNumberOfSequences();i++){
		//TODO : there was a pointer problem here
		scores.push_back(get_swps3_score_and_rc_cstyle(mat,&known_seqs->getSequence(i),&known_seqs->getSequence(i)));
	}

	for (int i=0;i<seqs.size();i++){
		if(i%reports == 0){
			cout << i << endl;
		}
		double maxide = 0;
		bool rc = false;
		for (int j=0;j<known_seqs->getNumberOfSequences();j++){
			bool trc = false;
			//TODO : there was a pointer problem here
			int ret = get_swps3_score_and_rc_cstyle(mat,&known_seqs->getSequence(j), & seqs[i]);
			double tsc = double(ret)/double(scores[j]);
			Sequence * rev = SequenceTools::reverse(seqs[i]);
			Sequence *comp = SequenceTools::complement(*rev);
			delete(rev);
			//TODO : there was a pointer problem here
			int retrc = get_swps3_score_and_rc_cstyle(mat,&known_seqs->getSequence(j), comp);
			delete(comp);
			if(retrc > ret){
				trc = true;
				tsc = double(retrc)/double(scores[j]);
			}
			if (tsc > maxide){
				maxide = tsc;
				rc = trc;
			}
		}//cout << maxide << endl;
		if (maxide >= identity){
			keep_seqs.push_back(seqs[i]);
			keep_rc.push_back(rc);
		}
	}
	my_data->keep_seqs = keep_seqs;
	my_data->keep_rc = keep_rc;
}


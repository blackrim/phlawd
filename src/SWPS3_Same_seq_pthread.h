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
 *  Same_seq_pthread.h
 */

#ifndef _SWPS3_SAME_SEQ_PTHREAD_H_
#define _SWPS3_SAME_SEQ_PTHREAD_H_

#include <string>
#include <vector>
#include <stdint.h>

using namespace std;

#include "sequence.h"
#include "DBSeq.h"
#include "utils.h"

struct SWPS3_thread_data{
	int thread_id;
	vector<DBSeq> seqs;
	vector<DBSeq> keep_seqs;
	vector<bool> keep_rc;
	int reports;
	double coverage;
	double identity;
	vector<Sequence> * known_seqs;
};

typedef int8_t * SBMatrix;
//took out the const
int get_swps3_score_and_rc_cstyle(SBMatrix mat, Sequence * inseq1, Sequence * inseq2);

void * SWPS3_Same_seq_pthread_go(void *threadarg);

#endif

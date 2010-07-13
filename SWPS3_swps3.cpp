/** \file swps3.c
 *
 * Main procedure and multi-threading code.
 */
/*
 * Copyright (c) 2007-2008 ETH Zürich, Institute of Computational Science
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SWPS3_swps3.h"
#include "SWPS3_matrix.h"
#include "SWPS3_fasta.h"
#include "SWPS3_DynProgr_scalar.h"
#ifdef SSE2
#include "SWPS3_DynProgr_sse_byte.h"
#include "SWPS3_DynProgr_sse_short.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/select.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <Seq/Sequence.h>
#include <Seq/containers>
#include <Seq/ioseq>
#include <Seq/alphabets>
using namespace bpp;

static ssize_t write_all(int fd, const void *data, size_t len) {
	size_t sent = 0;
	while(sent < len) {
		ssize_t res;
		res = write(fd,(const int8_t*)data+sent,len-sent);
		switch(res) {
		case 0:
			return sent;
		case -1:
			return -1;
		default:
			sent += res;
		}
	}
	return sent;
}

static ssize_t read_all(int fd, void *buf, size_t len) {
	size_t recv = 0;
	while(recv < len) {
		ssize_t res;
		res = read(fd,(int8_t*)buf+recv,len-recv);
		switch(res) {
		case 0:
			return recv;
		case -1:
			return -1;
		default:
			recv += res;
		}
	}
	return recv;
}


/*
 * switch( argv[i][1] ){
				case 'h':
					matrixFile = NULL;
					i = argc; break;
				case 's':
					type = SCALAR;
					break;
				case 't':
					options.threshold = atoi( argv[++i] );
					break;
				case 'i':
					options.gapOpen = atoi( argv[++i] );
					break;
				case 'e':
					options.gapExt = atoi( argv[++i] );
					break;
				case 'j':
					threads = atoi( argv[++i] );
					break;
				default:
					matrixFile = NULL;
					i = argc; break;
			}
 */
/*
 * have to send the matrix, which should be read in first
 * using SBMatrix swps3_readSBMatrix( char * filename );
 * long before sending to this
 *
 * send along the known seqs in the first and the query in
 * the second
 *
 * return the max score
 * query = known
 * db = test
 */
#include "SWPS3_fasta.h"
int swps3_maxscores ( SBMatrix matrix , const Sequence * known, const Sequence * test){
	int i,queryLen;
	const char * query;
	//SWType type = SSE2;
	Options options = {-12,-2,DBL_MAX};

	char * x1;
	x1 = (char*)malloc(sizeof(char)*known->size());
	memcpy(x1, known->toString().c_str(), known->size());
	swps3_translateSequence(x1,known->size(),NULL);
	query = x1;
	//query = known->toString().c_str();
	queryLen = known->size();
	double score = 0;
#ifdef SSE2
	ProfileByte  * profileByte = swps3_createProfileByteSSE( query, queryLen, matrix );
	ProfileShort * profileShort = swps3_createProfileShortSSE( query, queryLen, matrix );
#endif
	int dbLen;
	const char * db;

	char * x2;
	x2 = (char*)malloc(sizeof(char)*test->size());
	memcpy(x2, test->toString().c_str(), test->size());
	swps3_translateSequence(x2,test->size(),NULL);
	db = x2;
	//db=test->toString().c_str();
	dbLen = test->size();

#ifdef DEBUG
	for(i=0; i<queryLen; ++i) printf("\t%c",query[i]);
	printf("\n");
#endif

#ifdef SSE2
	//if(type == SSE2) {
		if( (score = swps3_alignmentByteSSE( profileByte, db, dbLen, &options )) >= DBL_MAX ) {
			score = swps3_alignmentShortSSE( profileShort, db, dbLen, &options );
			//assert(score >= 250 && "score too low");
		}
	//}
#else
	//if(type == SCALAR) {
		/*
		 * doesn't work!!
		 */
			double dmatrix[MATRIX_DIM*MATRIX_DIM];
			for(i=0;i<MATRIX_DIM*MATRIX_DIM;++i) dmatrix[i]=matrix[i];
			score = swps3_alignScalar( dmatrix, query, queryLen, db, dbLen, &options);
	//}
#endif


#ifdef SSE2
		swps3_freeProfileByteSSE( profileByte );
		swps3_freeProfileShortSSE( profileShort );
#endif
	free(x1);free(x2);
	return int(score);
}


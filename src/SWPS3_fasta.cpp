/** \file fasta.c
 *
 * Routines for reading FASTA databases.
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

#include "SWPS3_fasta.h"
#include "SWPS3_swps3.h"
#include "SWPS3_debug.h"
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif /* HAVE_MALLOC_H */
#include <string.h>
#include <errno.h>

EXPORT FastaLib * swps3_openLib( char * filename ){
	FastaLib * lib = NULL;
	FILE * fp;
	int len;
	if ( (fp = fopen( filename, "r" )) ){
#ifdef HAVE_MALLOC_H
		lib = memalign( 16, sizeof( FastaLib ) );
#else
		lib = (FastaLib*) malloc( sizeof( FastaLib ) );
#endif /* HAVE_MALLOC_H */
		lib->fp = fp;
		lib->name = "";
	}
	else{
		error("Fasta: %s\n",strerror(errno));
	}
	swps3_readNextSequence( lib, &len );
	if (len)
		rewind( lib->fp );
	return lib;
}
EXPORT char * swps3_readNextSequence( FastaLib * lib, int * len ){
	char * pres = lib->readBuffer;
	if ( feof( lib->fp ) )
		return NULL;
	if( (*pres = fgetc(lib->fp)) != '>' ) {
	   warning("Missing comment line, trying to continue anyway\n");
	   lib->name = "";
	   lib->data = pres++;
	   goto readseq;
	} else {
           fgets( pres, sizeof(lib->readBuffer)-(pres-lib->readBuffer), lib->fp );
	   lib->name = pres;
	   for (;*pres && *pres!='\n';pres++);
	   *pres++ = '\0';
	   while ((long)pres&0xf) *pres++ = '\0';
	   lib->data = pres;
	}
	while ( (*pres = fgetc( lib->fp) ) != '>' ){
readseq:
                if( fgets( pres+1, sizeof(lib->readBuffer)-(pres+1-lib->readBuffer), lib->fp ) == 0 ) goto finish;
	        for (;*pres && *pres!='\n';pres++)
			if ('A'>*pres || *pres>'Z'){
				error("Invalid character in input sequence '%c'\n", *pres);
			}
	}
	*pres = '\0';
	ungetc( '>', lib->fp );
finish:
	if (len)
		*len = pres - lib->data;
        bool isValid;
	isValid = swps3_translateSequence(lib->data,pres - lib->data,NULL);
	return lib->data;
}
EXPORT char * swps3_getSequenceName( FastaLib * lib ){
        return lib->name;
}
EXPORT void swps3_closeLib( FastaLib * lib ){
	free( lib );
}
EXPORT bool swps3_translateSequence(char * sequence, int seqLen, char table[256]) {
        /* edit sequence in place, replacing each char with its coding from the
         * lookup table if such a table is supplied. if not then convert to int
         * by subtracting 'A' from the char. Finally, validate the char. Return
         * true if all chars in the sequence are valid, false if not. */  
        int i;
        int lastBadChar = -1;
        for(i=0; i<seqLen && sequence[i]!='\n' && sequence[i]!='\0'; ++i) {
                // translate this char to integer code
                if(table) sequence[i] = table[(int)sequence[i]];
                else sequence[i] -= 'A';
                if(sequence[i] < 0 || sequence[i] >= MATRIX_DIM) {
                        // the integer value is outside the values expected
                        // (for capital chars this is 0-25. at the time of writing, MATRIX_DIM = 26)
                        lastBadChar = i;
                        printf("Invalid character in input sequence at position %d\n", lastBadChar);
                }
        }
        if (lastBadChar < 0) return true;
        else return false;
}

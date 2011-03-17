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
 * SmithWatermanGotoh.cpp
 */

#include <vector>
#include <string>
#include <limits>
#include <stdio.h>
#include <iostream>

using namespace std;

#include "SmithWatermanGotoh.h"
#include "utils.h"

const int EDNAFULL[][16] = {

		{5,-4,-4,-4,-4, 1, 1,-4,-4, 1,-4,-1,-1,-1,-2,-4},
		{-4, 5,-4,-4,-4, 1,-4, 1, 1,-4,-1,-4,-1,-1,-2, 5},
		{-4,-4, 5,-4, 1,-4, 1,-4, 1,-4,-1,-1,-4,-1,-2,-4},
		{-4,-4,-4, 5, 1,-4,-4, 1,-4, 1,-1,-1,-1,-4,-2,-4},
		{-4,-4, 1, 1,-1,-4,-2,-2,-2,-2,-1,-1,-3,-3,-1,-4},
		{ 1, 1,-4,-4,-4,-1,-2,-2,-2,-2,-3,-3,-1,-1,-1, 1},
		{ 1,-4, 1,-4,-2,-2,-1,-4,-2,-2,-3,-1,-3,-1,-1,-4},
		{-4, 1,-4, 1,-2,-2,-4,-1,-2,-2,-1,-3,-1,-3,-1, 1},
		{-4, 1, 1,-4,-2,-2,-2,-2,-1,-4,-1,-3,-3,-1,-1, 1},
		{ 1,-4,-4, 1,-2,-2,-2,-2,-4,-1,-3,-1,-1,-3,-1,-4},
		{-4,-1,-1,-1,-1,-3,-3,-1,-1,-3,-1,-2,-2,-2,-1,-1},
		{-1,-4,-1,-1,-1,-3,-1,-3,-3,-1,-2,-1,-2,-2,-1,-4},
		{-1,-1,-4,-1,-3,-1,-3,-1,-3,-1,-2,-2,-1,-2,-1,-1},
		{-1,-1,-1,-4,-3,-1,-1,-3,-1,-3,-2,-2,-2,-1,-1,-1},
		{-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2},
		{-4, 5,-4,-4,-4, 1,-4, 1, 1,-4,-1,-4,-1,-1,-2, 5}

};

const char IDENTITY	= '|';
const char SIMILARITY	= ':';
const char GAP		= ' ';
const char MISMATCH	= '.';
const char ALGAP = '-';
const unsigned int STOP = 0;
const unsigned int LEFT = 1;
const unsigned int DIAGONAL = 2;
const unsigned int UP = 3;
/**
 * Aligns two sequences by Smith-Waterman algorithm
 *
 * @param s1
 *            sequence #1
 * @param s2
 *            sequence #2
 * @param matrix
 *            scoring matrix
 * @param o
 *            open gap penalty
 * @param e
 *            extend gap penalty
 * @return alignment object contains the two aligned sequences, the
 *         alignment score and alignment statistics
 * @see Sequence
 * @see Matrix
 */


alignment * alignSWG(string s1, string s2, float o, float e) {
	int m = s1.length() + 1;
	int n = s2.length() + 1;

	vector<char> pointers(m * n);

	// Initializes the boundaries of the traceback matrix to STOP.
	for (int i = 0, k = 0; i < m; i++, k += n) {
		pointers[k] = STOP;
	}
	for (int j = 1; j < n; j++) {
		pointers[j] = STOP;
	}

	vector<short> sizesOfVerticalGaps(m * n);
	vector<short> sizesOfHorizontalGaps(m * n);
	for (int i = 0, k = 0; i < m; i++, k += n) {
		for (int j = 0; j < n; j++) {
			sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
		}
	}

	cell * ce = construct(s1, s2, o, e, &pointers,&sizesOfVerticalGaps, &sizesOfHorizontalGaps);
	alignment * align = traceback(s1, s2, &pointers, ce, &sizesOfVerticalGaps, &sizesOfHorizontalGaps);
	align->open = o;
	align->extend = e;
	delete ce;
	return align;
}

/**
 * Constructs directions matrix for the traceback
 *
 * @param s1
 *            sequence #1
 * @param s2
 *            sequence #2
 * @param matrix
 *            scoring matrix
 * @param o
 *            open gap penalty
 * @param e
 *            extend gap penalty
 * @return The cell where the traceback starts.
 */

cell * construct(string s1, string s2, float o,
		float e, vector<char> * pointers, vector<short> * sizesOfVerticalGaps,
		vector<short> * sizesOfHorizontalGaps) {

	int m = s1.length() + 1;
	int n = s2.length() + 1;

	float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
	vector<float> g(n); // score if xi aligns to a gap after yi
	float h; // score if yi aligns to a gap after xi
	vector<float> v(n); // best score of alignment x1...xi to y1...yi
	float vDiagonal = 0;

	g[0] = -std::numeric_limits<float>::infinity();
	h = -std::numeric_limits<float>::infinity();
	v[0] = 0;

	for (int j = 1; j < n; j++) {
		g[j] = -std::numeric_limits<float>::infinity();
		v[j] = 0;
	}

	float similarityScore, g1, g2, h1, h2;

	cell * ce = new cell();

	for (int i = 1, k = n; i < m; i++, k += n) {
		h = -std::numeric_limits<float>::infinity();
		vDiagonal = v[0];
		for (int j = 1, l = k + 1; j < n; j++, l++) {
			similarityScore = EDNAFULL[matrixTranslate(&s1[i - 1])][matrixTranslate(&s2[j - 1])];
			// Fill the matrices
			f = vDiagonal + similarityScore;

			g1 = g[j] - e;
			g2 = v[j] - o;
			if (g1 > g2) {
				g[j] = g1;
				sizesOfVerticalGaps->at(l) = (short) (sizesOfVerticalGaps->at(l - n) + 1);
			} else {
				g[j] = g2;
			}

			h1 = h - e;
			h2 = v[j - 1] - o;
			if (h1 > h2) {
				h = h1;
				sizesOfHorizontalGaps->at(l) = (short) (sizesOfHorizontalGaps->at(l - 1) + 1);
			} else {
				h = h2;
			}

			vDiagonal = v[j];
			v[j] = maximum(f, g[j], h, 0);

			// Determine the traceback direction
			if (v[j] == 0) {
				pointers->at(l) = STOP;
			} else if (v[j] == f) {
				pointers->at(l) = DIAGONAL;
			} else if (v[j] == g[j]) {
				pointers->at(l) = UP;
			} else {
				pointers->at(l) = LEFT;
			}

			// Set the traceback start at the current cell i, j and score
			if (v[j] > ce->score) {
				ce->row = i;
				ce->col = j;
				ce->score = v[j];
			}
		}
	}
	return ce;
}

/**
 * Returns the alignment of two sequences based on the passed array of
 * pointers
 *
 * @param s1
 *            sequence #1
 * @param s2
 *            sequence #2
 * @param m
 *            scoring matrix
 * @param cell
 *            The cell where the traceback starts.
 * @return {@link Alignment}with the two aligned sequences and alignment
 *         score.
 * @see Cell
 * @see Alignment
 */
alignment * traceback(string s1, string s2,
		vector<char> * pointers, cell * ce, vector<short> * sizesOfVerticalGaps,
		vector<short> * sizesOfHorizontalGaps) {

	int n = s2.length() + 1;

	alignment * align = new alignment();
	align->score = ce->score;

	int maxlen = s1.length() + s2.length(); // maximum length after the
											// aligned sequences

	string reversed1=""; // reversed sequence #1
	string reversed2=""; // reversed sequence #2
	string reversed3=""; // reversed markup

	int len1 = 0; // length of sequence #1 after alignment
	int len2 = 0; // length of sequence #2 after alignment
	int len3 = 0; // length of the markup line

	int identity = 0; // count of identitcal pairs
	int similarity = 0; // count of similar pairs
	int gaps = 0; // count of gaps

	char c1, c2;

	int i = ce->row; // traceback start row
	int j = ce->col; // traceback start col
	int k = i * n;

	bool stillGoing = true; // traceback flag: true -> continue & false
							   // -> stop

	while (stillGoing) {
		switch (pointers->at(k + j)) {
		case UP:
			for (int l = 0, len = sizesOfVerticalGaps->at(k + j); l < len; l++) {
				//reversed1[len1++] = s1[--i];
				//reversed2[len2++] = ALGAP;
				//reversed3[len3++] = GAP;
				reversed1.append(&s1[--i]);len1++;
				reversed2.append(&ALGAP);len2++;
				reversed3.append(&GAP);len3++;
				k -= n;
				gaps++;
			}
			break;
		case DIAGONAL:
			c1 = s1[--i];
			c2 = s2[--j];
			k -= n;
			//reversed1[len1++] = c1;
			//reversed2[len2++] = c2;
			reversed1.append(1,c1);len1++;
			reversed2.append(1,c2);len2++;
			if (c1 == c2) {
				//reversed3[len3++] = IDENTITY;
				reversed3.append(1,IDENTITY);len3++;
				identity++;
				similarity++;
			} else if (EDNAFULL[c1][c2] > 0) {
				//reversed3[len3++] = SIMILARITY;
				reversed3.append(1,SIMILARITY);len3++;
				similarity++;
			} else {
				//reversed3[len3++] = MISMATCH;
				reversed3.append(1,MISMATCH);len3++;
			}
			break;
		case LEFT:
			for (int l = 0, len = sizesOfHorizontalGaps->at(k + j); l < len; l++) {
				//reversed1[len1++] = ALGAP;
				//reversed2[len2++] = s2[--j];
				//reversed3[len3++] = GAP;
				reversed1.append(1,ALGAP);len1++;
				reversed2.append(1,s2[--j]);len2++;
				reversed3.append(1,GAP);len3++;
				gaps++;
			}
			break;
		case STOP:
			stillGoing = false;
		}
	}

	align->seq1 = (reverse(reversed1, len1));
	align->start1 = i;
	align->seq2 = (reverse(reversed2, len2));
	align->start2 = j;
	align->markup = (reverse(reversed3, len3));
	align->identity = (identity);
	align->gaps = (gaps);
	align->similarity = (similarity);
	cout << "+ " << gaps << endl;
	return align;
}

/**
 * Returns the maximum of 4 float numbers.
 *
 * @param a
 *            float #1
 * @param b
 *            float #2
 * @param c
 *            float #3
 * @param d
 *            float #4
 * @return The maximum of a, b, c and d.
 */
float maximum(float a, float b, float c, float d) {
	if (a > b) {
		if (a > c) {
			return a > d ? a : d;
		} else {
			return c > d ? c : d;
		}
	} else if (b > c) {
		return b > d ? b : d;
	} else {
		return c > d ? c : d;
	}
}

/**
 * Reverses an array of chars
 *
 * @param a
 * @param len
 * @return the input array of char reserved
 */
string reverse(string a, int len) {
	string b(a);
	for (int i = len - 1, j = 0; i >= 0; i--, j++) {
		b[j] = a[i];
	}
	return b;
}

int matrixTranslate(string ch){
	if ( "A")
		return 0;
	else if ( "T")
		return 1;
	else if ( "G")
		return 2;
	else if ( "C")
		return 3;
	else if ( "S")
		return 4;
	else if ( "W")
		return 5;
	else if ( "R")
		return 6;
	else if ( "Y")
		return 7;
	else if ( "K")
		return 8;
	else if ( "M")
		return 9;
	else if ( "B")
		return 10;
	else if ( "V")
		return 11;
	else if ( "H")
		return 12;
	else if ( "D")
		return 13;
	else if ( "N")
		return 14;
	else if ( "U")
		return 15;
}

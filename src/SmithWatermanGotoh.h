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

#ifndef _SMITHWATERMANGOTOH_H_
#define _SMITHWATERMANGOTOH_H_

#include <string>
#include <vector>

using namespace std;


struct cell{
int row;
int col;
float score;
};

struct alignment{
float score;
float open;
float extend;
string seq1;
string seq2;
int start1;
int start2;
string markup;
int identity;
int gaps;
int similarity;
};


alignment * alignSWG(string s1, string s2, float o, float e);
cell * construct(string s1, string s2, float o,
		float e, vector<char> * pointers, vector<short> * sizesOfVerticalGaps,
		vector<short> * sizesOfHorizontalGaps);
alignment * traceback(string s1, string s2,
		vector<char> * pointers, cell * ce, vector<short> * sizesOfVerticalGaps,
		vector<short> * sizesOfHorizontalGaps);
float maximum(float a, float b, float c, float d);
string reverse(string a, int len);
int matrixTranslate(string ch);
#endif

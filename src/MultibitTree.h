// MultibitTree.h
//
// Copyright (c) 2015
// Universitaet Duisburg-Essen
// Campus Duisburg
// Institut fuer Soziologie
// Prof. Dr. Rainer Schnell
// Lotharstr. 65
// 47057 Duisburg 
//
// This file is part of the R-Package "multibitTree".
//
// "multibitTree" is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// "multibitTree" is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with "multibitTree". If not, see <http://www.gnu.org/licenses/>.

#ifndef MULTIBITTREE_H
#define MULTIBITTREE_H

#include "Fingerprint.h"
#include "QueryResult.h"

// Objects of class MultibitTree contain a tree data structure for
// performing a fast Tanimoto search. The underlying array of prints
// will be sorted according to the tree structure.
//
// The MultibitTree data structure is based on the MultibitTree described in
// http://www.almob.org/content/5/1/9

typedef unsigned short ushort;

class MultibitTree {
	private:
	
	int mCardinality;		// common cardinality for this tree
	int mLeafLimit;			// the maximum number of fingerprints for which
					// no further sub-tree shall be calculated
	int mNBits;			// maximal size of Fingerprints in bits
	long long mLeafStart;		// start of MultibitTree in the array of Fingerprints
	long long mSize;		// length of MultibitTree in the array of Fingerprints
	long long mTreeSize;		// size of tree data structure
	long long mNodes;		// count of tree nodes
	ushort **mMatchBits;		// match-bit-list for each tree node
	ushort *mMatchBitsSize;		// total size of match bits for each tree node
	ushort *mMatchBitsZerosSize;	// size of zero match bits for each tree node
	long long *mLeftChild;		// left subtree node for each inner node
					// for leaf nodes this is reused for index into leaves list (start)
	long long *mRightChild;		// right subtree node for each inner node
					// for leaf nodes this is reused for index into leaves list (end)
	Fingerprint **mLeaves;		// pointer to an array of Fingerprints
	
	long long *mCounts;		// counter for match bit computation
	ushort *mMatchListZeros;	// temporary list for match bit computation
	ushort *mMatchListOnes;		// temporary list for match bit computation
	
	long long mCntXOR;		// statistic counter before XOR-check
	long long mCntTanimoto;		// statistic counter before tanimoto-check
	
	// recursively build a subtrees
	void buildNode(Fingerprint *usedBits, long long leafStart, long long leafEnd);

	// sort sub tree prints by best match bit
	long long splitLeavesHalf(long long leafStart, long long leafEnd);
	
	// recursively search sub tree
	void internalSearch (QueryResult *result, Fingerprint *queryPrint, long long node, int commonXOR, int AB, int queryUnmatched, int treeUnmatchedi, float minTanimoto);

	public:
	
	// constructor
	// create a new MultibitTree from an array of Fingerprints
	MultibitTree(Fingerprint **prints, long long leafStart, long long leafEnd, int nBits, int cardinality, int leafLimit);
	
	// destructor
	~MultibitTree();

	// perform a search for <query> that has <cardinality> filtered by <minTanimoto>
	// and add the result to <result>
	inline void search(QueryResult *result, Fingerprint *queryPrint, int cardinality, float minTanimoto) {
		internalSearch(result, queryPrint, 0, 0, cardinality + mCardinality, cardinality, mCardinality, minTanimoto);
	}
	
	// return tree size
	inline long long getSize() {
		return mSize;
	}
	
	// return XOR counter
	inline long long getCntXOR() {
		return mCntXOR;
	}
	
	// return tanimoto counter
	inline long long getCntTanimoto() {
		return mCntTanimoto;
	}

	// initialize XOR counter
	inline void initCntXOR() {
		mCntXOR = 0;
	}
	
	// inittialize tanimoto counter
	inline void initCntTanimoto() {
		mCntTanimoto = 0;
	}
};
#endif

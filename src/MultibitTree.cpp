// MultibitTree.cpp
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

#include <stdio.h>
#include <assert.h>

#include "MultibitTree.h"
#include "Misc.h"

#define LEAF_BIT 0x8000
#define BIT_MASK 0x7fff

// constructor
// create a new MultibitTree from an array of Fingerprints
// prints		pointer to array of Fingerprints
// leafStart		cluster starting position in prints
// leafEnd		end of cluster
// nBits		maximal size of Fingerprint in bits
// cardinality		cluster cardinality
// leafLimit		leaf limit for MultibitTree creation
MultibitTree::MultibitTree(Fingerprint **prints, long long leafStart, long long leafEnd, int nBits, int cardinality, int leafLimit) {
	Fingerprint *usedBits;
	
	// initialize member fields and allocate memory
	mSize = leafEnd - leafStart;
	mLeafStart = leafStart;
	mNBits = nBits;
	mCardinality = cardinality;
	mLeafLimit = leafLimit;
	mTreeSize = MAX(2 * mSize - 1, 1);		// max size of binary-tree with "size" leaves
	mMatchBits = new ushort*[mTreeSize];
	mMatchBitsSize = new ushort[mTreeSize];
	mMatchBitsZerosSize = new ushort[mTreeSize];
	mLeftChild = new long long[mTreeSize];
	mRightChild = new long long[mTreeSize];
	mLeaves = prints;

	mCounts = new long long[nBits];
	mMatchListOnes = new ushort[nBits];
	mMatchListZeros = new ushort[nBits];
	usedBits = new Fingerprint(nBits);
	mNodes = 0;
	mCntXOR = 0;
	mCntTanimoto = 0;

	// build tree
	buildNode(usedBits, leafStart, leafEnd);
	
	// free temporary data structures
	delete usedBits;
	delete[] mMatchListZeros;
	delete[] mMatchListOnes;
	delete[] mCounts;
}

// destructor
// delete all used resources
MultibitTree::~MultibitTree() {
	for (long long i = 0; i < mNodes; i++) {
		if ((mMatchBitsSize[i] & BIT_MASK) != 0) {
			delete[] mMatchBits[i];
		}
	}
	
	delete[] mMatchBits;
	delete[] mMatchBitsSize;
	delete[] mLeftChild;
	delete[] mRightChild;
}

// recursively build MultibitTree nodes and sub-nodes
void MultibitTree::buildNode(Fingerprint *usedBits, long long leafStart, long long leafEnd) {
	ushort listCountOnes;		// match bit counter for 1-bits
	ushort listCountZeros;		// match bit counter for 0-bits
	ushort listCount;		// match bit counter total; 
	long long nLeaves;		// size of sub-tree range in leaves
	long long middle;		// index for sub-tree splitting
	long long thisNode;		// index of current node
	ushort *list;			// temporary match bit buffer
	Fingerprint *clone;		// used bits clone for left sub-tree
	
	nLeaves = leafEnd - leafStart;

	thisNode = mNodes;
	mNodes++;
	
	assert(thisNode < mTreeSize);
	
	// count bits and find match-bits
	listCountOnes = 0;
	listCountZeros = 0;

	for (int i = 0; i < mNBits; i++) {
		mCounts[i] = 0;
		if (!usedBits->getBit(i)) {
			// for each unused bit count 1-bits within sub-tree range
			for (int leaf = leafStart; leaf < leafEnd; leaf++) {
				if (mLeaves[leaf]->getBit(i)) {
					mCounts[i]++;
				}
			}

			if (mCounts[i] == 0) {
				// if all bits are "0" store match-bit in zero-list
				mMatchListZeros[listCountZeros] = i;
				listCountZeros++;
				usedBits->setBit(i);
			} else if (mCounts[i] == nLeaves) {
				// if all bits are "1" store match-bit in one-list
				mMatchListOnes[listCountOnes] = i;
				listCountOnes++;
				usedBits->setBit(i);
			}
		}
	}

	// copy match-bits in a well-sized persistent array

	listCount = listCountOnes + listCountZeros;

	if (listCount > 0) {
		list = new ushort[listCount];
		for (int i = 0; i < listCountZeros; i++) {
			list[i] = mMatchListZeros[i];
		}
		for (int i = 0; i < listCountOnes; i++) {
			list[listCountZeros + i] = mMatchListOnes[i];
		}
		mMatchBits[thisNode] = list;
	}

	mMatchBitsSize[thisNode] = listCount;
	mMatchBitsZerosSize[thisNode] = listCountZeros;

	// TDB: this could be checked earlier, if the deletion below could be omitted
	if (nLeaves < mLeafLimit) {
		middle = -1; // leaf
	} else {
		// TDB: this could be done earlier, if the deletion below could be omitted
		middle = splitLeavesHalf(leafStart, leafEnd);
	}
	
	if (middle == -1) {
		// create leaf
		mMatchBitsSize[thisNode] = LEAF_BIT;
		// TBD: original implementation requires a deletion of match-bits for leaf-nodes
		// check, if using them here may be more reasonable, too
		delete[] mMatchBits[thisNode];
		mLeftChild[thisNode] = leafStart;
		mRightChild[thisNode] = leafEnd;
	} else {
		// create inner node

		// copy to temporary clone of used bits for left sub-tree
		clone = new Fingerprint(usedBits);

		// build left sub-tree
		mLeftChild[thisNode] = mNodes;
		buildNode(clone, leafStart, middle);
		
		// delete clone
		delete clone;
		
		// build right sub-tree
		mRightChild[thisNode] = mNodes;
		buildNode(usedBits, middle, leafEnd);
	}
}

// Clustering Strategy
// split leaves (clustering strategy: split-half)
long long MultibitTree::splitLeavesHalf(long long leafStart, long long leafEnd) {
	int bestBit;		// best bit for sub-tree clustering
	int bestDist;		// best (lowest) distance between 1-bits and half of the clusters size
	int ones;		// number of 1-bits in the cluster
	int currentDist;	// current distance between 1-bits and half of the clusters size
	long long left, right;	// left and right sort index for 0/1-sorting
	long long half, nLeaves;	// half and total size of leaf-cluster
	Fingerprint *leaf;	// temporary pointer for swapping Fingerprints

	// find best bit for sub-tree clustering
	bestBit = -1;
	bestDist = mNBits;
	nLeaves = leafEnd - leafStart;
	half = nLeaves / 2;

	// computes best bit by lowest distance between 1-bits and half of the clusters size
	for (int i = 0; i < mNBits; i++) {
		ones = mCounts[i];
		currentDist = ABS(ones - half);
		
		if ((ones != 0) && (ones != nLeaves) && (currentDist < bestDist)) {
			bestBit = i;
			bestDist = currentDist;
		}
	}

	if (bestBit == -1) {
		return -1; // no split => leaf
	}
		
	// sort and split sub-tree by bit
	// all Fingerprints with 0 at the bestBit-position are sorted to the left and
	// all Fingerprints with 1 at the bestBit-position are sorted to the right
	left = leafStart;
	right = leafEnd - 1;
	
	while (left < right) {
		if (!mLeaves[left]->getBit(bestBit)) {
			left++;
			continue;
		}

		if (mLeaves[right]->getBit(bestBit)) {
			right--;
			continue;
		}
		
		// swap leaves
		leaf = mLeaves[left];
		mLeaves[left] = mLeaves[right];
		mLeaves[right] = leaf;
	}

	return right;
}
	
// searching
// traverse the tree and visit only those sub-trees that don't surely underrun the tanimoto filter
void MultibitTree::internalSearch(QueryResult *result, Fingerprint *queryPrint, long long node, int commonXOR, int AB, int queryUnmatched, int treeUnmatched, float minTanimoto) {	
	int size;
	ushort *matchBitIdx;

	size = mMatchBitsSize[node];

	if (size & LEAF_BIT) {
		// if this is a leaf-node check each leaf-Fingerprint's tanimoto coefficient
		for (long long i = mLeftChild[node]; i < mRightChild[node]; i++) {
			Fingerprint *leaf = mLeaves[i];

			// increase statistic counter for XOR-hash estimation
			mCntXOR++;

			// check XOR-hash estimation
			if (queryPrint->tanimotoXOR(leaf, AB) >= minTanimoto) {
				// increase statistic counter for tanimoto calculation
				mCntTanimoto++;

				// check exact tanimoto condition and add to results if matches
				float tanimoto = queryPrint->tanimoto(leaf);

				// check exact tanimoto condition
				if (tanimoto >= minTanimoto) {
					// add matching leaf to QueryResult
					result->add(queryPrint->getId(), leaf, tanimoto);
				}
			}
		}
	} else {
		// if this is an inner node estimate minimal tanimoto-coefficient
		// by evaluating the match-bits
		
		int countOnes = 0;
		int countZeros = 0;

		int sizeZeros =  mMatchBitsZerosSize[node];
		matchBitIdx = mMatchBits[node];

		// count differences: count 1-bits in query for all 0-bits in match-bits
		for (int i = 0; i < sizeZeros; i++) {
			countOnes += queryPrint->getBit(matchBitIdx[i]);
		}

		// count differences: count 0-bits in query for all 1-bits in match-bits
		for (int i = sizeZeros; i < size; i++) {
			countZeros += queryPrint->getBit(matchBitIdx[i]) ^ 1;
		}

		commonXOR += countZeros + countOnes;	// total differences
		queryUnmatched -= countOnes;		// subtract 1-bits of query-bits covered by match-bits
		treeUnmatched -= countZeros;		// subtract 1-bits of tree-bits  covered by match-bits

		// compute and compare minimal tanimoto-coefficient for this sub-tree
		if (((float) MIN(queryUnmatched, treeUnmatched)) / (commonXOR + MAX(queryUnmatched, treeUnmatched)) >= minTanimoto) {
			// analyse sub-trees
			internalSearch(result, queryPrint, mLeftChild[node], commonXOR, AB, queryUnmatched, treeUnmatched, minTanimoto);
			internalSearch(result, queryPrint, mRightChild[node], commonXOR, AB, queryUnmatched, treeUnmatched, minTanimoto);
		}
	}
}

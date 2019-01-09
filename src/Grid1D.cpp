// Grid1D.cpp
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

#include "Grid1D.h"

// constructor:
//
// sort Fingerprints by cardinality and
// create a MultibitTree for each cardinality
//
// prints	: pointer on Fingerprint array
// size		: size of <prints>
// nBits	: maximal size of Fingerprints in <prints>
// threads	: number of parallel threads passed to ThreadPool
// leafLimit	: leaf limit parameter passed to all MultibitTrees

Grid1D::Grid1D(Fingerprint **prints, long long size, int nBits, int threads, int leafLimit) {
	int count[nBits + 1];		// cardnality cluster counter
	int pos[nBits + 1];		// destination positions for each cluster
	Fingerprint *swap1, *swap2;	// helper pointer for sorting
	int card;			// helper variable for current cardinality

	mWorkerPool = new ThreadPool(threads);
	
	mNBits = nBits;
	mSize = size;
	mSizeLastSearch = 0;
	mBuckets = new MultibitTree*[nBits + 1];

	// sort prints by cardinality
	// this can be done in linear time because we have a limited number of clusters

	// initialize counters
	for (int i = 0; i < (nBits + 1); i++) {
		count[i] = 0;
	}

	// run through all prints and count the occurence of all possible cardinalities
	for (long long i = 0; i < size; i++) {
		count[prints[i]->cardinality()]++;
	}
	
	// pre-calculate the sorting-destination position for each cardinality
	// re-use count-array as cluster destination index for each cardinality
	pos[0] = 0;
	for (int i = 1; i < (nBits + 1); i++) {
		pos[i] = pos[i - 1] + count[i - 1];
		count[i - 1] = pos[i - 1];
	}
	count[nBits] = pos[nBits];
	
	// swap each print that is not in the right cardinality cluster
	for (long long i = 0; i < (size - 1); i++) {
		swap1 = prints[i];
		card = swap1->cardinality();

		// do swapping only if print is not in the right cluster
		// TODO: (i >= count[card]) may be redundant?
		
		if ((i < pos[card]) || (i >= count[card])) {
			// re-swap until the current destination equals i
			while (i != count[card]) {
				// store old print from destination
				swap2 = prints[count[card]];
				// store new print to its destination
				prints[count[card]] = swap1;
				// increment cluster destination index
				count[card]++;
				// use old print as new candidate for swapping
				swap1 = swap2;
				card = swap1->cardinality();
			}
			prints[i] = swap1;
			count[card]++;
		}
	}

	// create a MultibitTree for each cardinality cluster
	for (int i = 0; i < (nBits + 1); i++) {
		if (pos[i] < count[i]) {
			mWorkerPool->createMultibitTree(&mBuckets[i], prints, pos[i], count[i], nBits, i, leafLimit);
		} else {
			mBuckets[i] = NULL;
		}
	}

	// wait for running threads
	mWorkerPool->wait();
}

// destructor
// delete all used resources

Grid1D::~Grid1D() {
	for (int i = 0; i < mNBits; i++) {
		if (mBuckets[i]) {
			delete mBuckets[i];
		}
	}

	delete[] mBuckets;
	delete mWorkerPool;
}

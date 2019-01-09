// Grid1D.h
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

#ifndef GRID1D_H
#define GRID1D_H

#include <math.h>
#include "Fingerprint.h"
#include "MultibitTree.h"
#include "ThreadPool.h"

// Objects of class Grid1D hold an array of instances of the class MultibitTree.
// In each MultibitTree all Fingerprints of the same cardinality are stored.
// Knowing the queries cardinality and the Tanimoto coefficient a search can
// be reduced on a relevant range of MultibitTrees.
// The Grid1D also uses the class ThreadPool for concurrently work on different
// MultibitTrees.
//
// The Grid1D data structure is based on the kDGrid described in
// http://www.almob.org/content/5/1/9

class Grid1D {
	private:

	MultibitTree **mBuckets;	// array of MultibitTrees
	int mNBits;			// maximal size of Fingerprints
	long long mSize;		// size of Fingerprint-array used by mBuckets
	long long mSizeLastSearch;	// for statistics
	ThreadPool *mWorkerPool;	// ThreadPool for concurrency
	
	public:
	
	// constructor
	Grid1D(Fingerprint **prints, long long size, int nBits, int threads, int leafLimit);

	// destructor	
	~Grid1D();

	// perform a search for <query> and <minTanimoto> and add the result to <result>
	// parallelise by buckets
	inline void search(QueryResult *result, Fingerprint *query, float minTanimoto) {
		int min, max, card;

		card = query->cardinality();
		
		// search only in MultibitTrees with suitable cardinality
		min = (int) ceil((minTanimoto * card));
		max = MIN((int) (1.0 / minTanimoto * card) + 1, mNBits + 1);
		
		for (int i = min; i < max; i++) {
			if (mBuckets[i]) {
				mWorkerPool->searchMultibitTree(mBuckets[i], result, query, card, minTanimoto);
			}
		}
		
		// wait for running threads
		mWorkerPool->wait();
	}

	// perform a search for <query and <minTanimoto> and add the result to <result>
	// start search as one thread and return
	inline void searchAsync(QueryResult *result, Fingerprint *query, float minTanimoto) {
		int min, max, card;

		card = query->cardinality();
		
		// search only in MultibitTrees with suitable cardinality
		min = (int) ceil((minTanimoto * card));
		max = MIN((int) (1.0 / minTanimoto * card) + 1, mNBits + 1);
				
		mWorkerPool->searchMultibitTreeRange(mBuckets, min, max, result, query, card, minTanimoto);
	}

	// wait for running threads();
	inline void wait() {
		mWorkerPool->wait();
	}

	// init Statistic Values;
	inline void initStatistics() {
		for (int i = 0; i <= mNBits; i++) {
			if (mBuckets[i]) {
				mBuckets[i]->initCntXOR();
				mBuckets[i]->initCntTanimoto();
      			}
    		}
	}

	// get Statistics of last search
	inline void getStatistics(double *valuesPtr, double *percentsPtr) {
		long long cntX = 0;
		long long cntT = 0;

		for (int i = 0; i <= mNBits; i++) {
			if (mBuckets[i] != NULL) {
			  cntX += mBuckets[i]->getCntXOR(); 				
			  cntT += mBuckets[i]->getCntTanimoto();
			} 				
		}

		valuesPtr[0] = (double)cntX;
		valuesPtr[1] = (double)cntT;
		valuesPtr[2] = (double)(mSize * mSizeLastSearch);
    
		percentsPtr[0] = (double)cntX / (mSize * mSizeLastSearch) * 100;
		percentsPtr[1] = (double)cntT / (mSize * mSizeLastSearch) * 100;
		percentsPtr[2] = 100.0;
	}
	
	// set size of last search;
	inline void setSizeLastSearch(long long sls) {
		mSizeLastSearch = sls;
	}

	// get size of last search;
	inline long long getSizeLastSearch() {
		return mSizeLastSearch;
	}

};
#endif

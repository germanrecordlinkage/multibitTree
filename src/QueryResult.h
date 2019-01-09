// QueryResult.h
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

#ifndef QUERYRESULT_H
#define QUERYRESULT_H

#include <pthread.h>
#include "Fingerprint.h"

// Instances of QueryResultNode store one Fingerprint and
// the corresponding Tanimoto coefficient within a sorted
// tree. If sorting is disabled, QueryResultNodes are simply
// used as nodes of a linked list.
 
typedef struct QueryResultNodeStruct {
	char *mQueryId;				// ID of query Fingerprint
	Fingerprint *mPrint;			// pointer to matching Fingerprint
	float mTanimoto;			// corresponding Tanimoto coefficient
	struct QueryResultNodeStruct *mLeft;	// left subtree (higher Tanimoto coeffs)
	struct QueryResultNodeStruct *mRight;	// right subtree (lower Tanimoto coeffs)
} QueryResultNode;

// Objects of class QueryResult store search results,
// which consist of a Fingerprint and the computed
// Tanimoto coefficient. The QueryResult is designed as
// a sorted tree. On construction the QueryResult can be
// flagged to sort the results while collecting them by
// the Tanimoto coefficient or not. 
// if a result file is specified the search results
// will be stored directly into this file without
// using internal memory

class QueryResult {
	private:
	
	QueryResultNode *mRootNode;		// root node of sorted tree or linked list
	int mSize;				// size of result set
	int mSort;				// flag if the query results have ro be sorted
	pthread_mutex_t mAddMutex;		// thread save insertion mutex
	FILE *mResultFile;			// file descriptor for the optional result file
	const char *mSeperator;			// column seperator for csv output

	void deleteNode(QueryResultNode *node);	// delete a node and all its subnodes

	public:

	// constructor
	QueryResult(int sort, FILE *resultFile, const char *seperator);

	// destructor
	~QueryResult();
	
	// add a Fingerprint and the corresponding Tanimoto coefficient to the query result
	void add(char *queryId, Fingerprint *print, float tanimoto);
	
	// return root node for reading
	inline QueryResultNode *getRootNode() {
		return mRootNode;
	}

	// return size of result set
	inline int getSize() {
		return mSize;
	}
};
#endif

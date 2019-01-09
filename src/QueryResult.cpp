// QueryResult.cpp
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

#include "QueryResult.h"

// constructor
QueryResult::QueryResult(int sort, FILE *resultFile, const char *seperator) {
	mSize = 0;
	mRootNode = NULL;
	mSort = sort;
	mResultFile = resultFile;
	mSeperator = seperator;

	pthread_mutex_init(&mAddMutex, NULL);
}

// destructor
QueryResult::~QueryResult() {
	QueryResultNode *node, *nextNode;
	if (mSort) {
		deleteNode(mRootNode);		// discard the whole result set recursivly
	} else {
		// discard the whole result iteratively for large scale results
		node = mRootNode;
		while (node != NULL) {
			if (node->mQueryId != NULL) {
				delete[] node->mQueryId;
			}
			nextNode = node->mRight;
			delete node;
			node = nextNode;
		}
	}
}

// delete a node and all its subnodes recursively
void QueryResult::deleteNode(QueryResultNode *node) {
	if (node != NULL) {
		if (node->mQueryId != NULL) {
			delete[] node->mQueryId;
		}
		deleteNode(node->mLeft);
		deleteNode(node->mRight);
		delete node;
	}
}

// add a Fingerprint and the corresponding Tanimoto coefficient to the query result
void QueryResult::add(char *queryId, Fingerprint *print, float tanimoto) {
	// lock mutex
	pthread_mutex_lock(&mAddMutex);	

	// check if a results file is specified
	if (mResultFile == NULL) {
		// if not, collect results in memory
		QueryResultNode *newNode;
		QueryResultNode **nextNodePtr;
	
		// create new node
		newNode = new QueryResultNode;

		// copy query id, print pointer and tanimoto
		if (queryId != NULL) {
			newNode->mQueryId = new char[strlen(queryId) + 1];
			strcpy(newNode->mQueryId, queryId);
		} else {
			newNode->mQueryId = NULL;
		}

		newNode->mPrint = print;
		newNode->mTanimoto = tanimoto;
		newNode->mLeft = NULL;
		newNode->mRight = NULL;

		if (mSort) {
			// if query results have to be sorted, find the right location
			// to insert the new node
			nextNodePtr = &mRootNode;

			while (*nextNodePtr != NULL) {
				if ((*nextNodePtr)->mTanimoto >= tanimoto) {
					nextNodePtr = &((*nextNodePtr)->mRight);
				} else {
					nextNodePtr = &((*nextNodePtr)->mLeft);
				}
			}

			// insert new node
			*nextNodePtr = newNode;
		} else {
			// if query results need not to be sorted
			// insert at end of list
			if (mRootNode == NULL) {
				mRootNode = newNode;
			} else {
				mRootNode->mLeft->mRight = newNode;
			}

			mRootNode->mLeft = newNode; // use mLeft as pointer to end node
		}
	} else {
		// write results to file
		fprintf(mResultFile, "%s%s%s%s%.7f\n", queryId, mSeperator, print->getId(), mSeperator, tanimoto);
	}
	
	mSize++;

	// unlock mutex
	pthread_mutex_unlock(&mAddMutex);
}

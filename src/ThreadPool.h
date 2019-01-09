// ThreadPool.h
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

#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <pthread.h>
#include "QueryResult.h"
#include "Fingerprint.h"
#include "MultibitTree.h"

// forward declaration
class ThreadPool;

// Instances of threadDataType store task information
// for a running thread. The main thread dispatches tasks
// by filling in the task information into the task's slot
// and signaling the condition to start.
typedef struct threadDataStruct {
	int slot;			// the thread's slot
        ThreadPool *pool;		// pointer on responsible ThreadPool
	int task;			// task
	pthread_mutex_t mutex;		// mutex for signal handling
	pthread_cond_t condition;	// condition for signal handling
} threadDataType;

// Instances of createArgumentsType hold the parameters
// for performing the creation of a MultibitTree.
typedef struct createArgumentsStruct {
        MultibitTree **tree;		// address where to store the new MultibitTree pointer
        Fingerprint **prints;		// pointer to array of Fingerprints
        int leafStart;			// cluster starting position in prints
        int leafEnd;			// end of cluster
        int nBits;			// maximal size of Fingerprint in bits
        int cardinality;		// cluster cardinality
        int leafLimit;			// leaf limit for MultibitTree creation
} createArgumentsType;

// Instances of searchArgumentsType hold the parameters
// for searching in a MultibitTree.
typedef struct searchArgumentsStruct {
        MultibitTree *tree;		// pointer to the MultibitTree to search
        QueryResult *result;		// QueryResult for storing the results
        Fingerprint *query;		// query Fingerprint to search for
        int cardinality;		// cardinality of query
        float minTanimoto;		// filter criteria
} searchArgumentsType;

// Instances of searchRangeArgumentsType hold the parameters
// for searching in an array of MultibitTrees.
typedef struct searchRangeArgumentsStruct {
        MultibitTree **trees;		// pointer to the MultibitTree array to search
        int min;			// start of the range in the array
        int max;			// end of the range in the array
        QueryResult *result;		// QueryResult for storing the results
        Fingerprint *query;		// query Fingerprint to search for
        int cardinality;		// cardinality of query
        float minTanimoto;		// filter criteria
} searchRangeArgumentsType;

// Instances of ThreadPool hold a set threads that can concurrently
// perform task. ThreadPool dispatches a new task to the next free
// thread and returns. If all threads are working, ThreadPool waits
// until the task can be dispatched.
class ThreadPool {
	private:

	pthread_t *mPool;				// array of threads
	threadDataType *mThreadData;			// array of task information per thread
	createArgumentsType *mCreateArgs;		// array of arguments for task "create"
	searchArgumentsType *mSearchArgs;		// array of arguments for task "search"
	searchRangeArgumentsType *mSearchRangeArgs;	// array of arguments for task "searchRange"
	
	pthread_mutex_t mSlotStackMutex;	// mutex for signal handling
	pthread_cond_t mSlotStackCondition;	// condition for signal handling

	int mPoolSize;				// number of managed threads
	int *mSlotStack;			// stack of free threads
	int mSlotStackPosition;			// stack position

	int getSlot();				// lock next free thread
	void startSlot(int task, int slot);	// start locked thread to perform task
	void releaseSlot(int slot);		// unlock thread after completing a task

	public:

	// constructor
	ThreadPool(int size);			// create a new ThreadPool with size threads

	// destructor
	~ThreadPool();				// stop all threads and discard allocated resources
	
	void worker(int slot);			// thread main loop for retrieving and performing tasks

	// dispatch a task to create a new MultibitTree
	void createMultibitTree(MultibitTree **tree, Fingerprint **prints, int leafStart, int leafEnd, int nBits, int cardinality, int leafLimit);
	
	// dispatch a task to search in a MultibitTree
	void searchMultibitTree(MultibitTree *tree, QueryResult *result, Fingerprint *query, int cardinality, float minTanimoto);

	// dispatch a task to search in an array of MultibitTrees
	void searchMultibitTreeRange(MultibitTree **trees, int min, int max, QueryResult *result, Fingerprint *query, int cardinality, float minTanimoto);

	// wait until all threads have completed
	void wait();
};
#endif

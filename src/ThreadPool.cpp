// ThreadPool.cpp
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

#include <assert.h>
#include "ThreadPool.h"

// thread wrapper that is compatible to pthread-API and calls the
// worker member function of ThreadPool after thread creation
void *threadWrapper(void *p) {
	threadDataType *threadData;
	
	threadData = ((threadDataType*)p);
	threadData->pool->worker(threadData->slot);
	
	return NULL;
}

// thread main loop for retrieving and performing tasks
void ThreadPool::worker(int slot) {
	int running = true;
	int *task = &(mThreadData[slot].task);
	pthread_mutex_t *mutex = &(mThreadData[slot].mutex);
	pthread_cond_t *condition = &(mThreadData[slot].condition);
	
	while (running) {
		// wait for task

		pthread_mutex_lock(mutex);

		while (*task == 0) {
			pthread_cond_wait(condition, mutex);
		}
		
		pthread_mutex_unlock(mutex);
		
		// perform task

		if (*task == 1) {
			// create a new MultibitTree
			createArgumentsType *args = &(mCreateArgs[slot]);
			*(args->tree) = new MultibitTree(args->prints, args->leafStart, args->leafEnd, args->nBits, args->cardinality, args->leafLimit);
		} else if (*task == 2) {
			// search in a MultibitTree
			searchArgumentsType *args = &(mSearchArgs[slot]);
			args->tree->search(args->result, args->query, args->cardinality, args->minTanimoto);
		} else if (*task == 4) {
			// search in an array of MultibitTrees
			searchRangeArgumentsType *args = &(mSearchRangeArgs[slot]);
			for (int i = args->min; i < args->max; i++) {
				if (args->trees[i]) {
					args->trees[i]->search(args->result, args->query, args->cardinality, args->minTanimoto);
				}
			}
		} else if (*task == 3) {
			// stop thread
			running = false;
		}

		// signal the thread to be ready for a new task
		releaseSlot(slot);
	}
}

// constructor
// create a new ThreadPool with size threads
ThreadPool::ThreadPool(int size) {
	pthread_attr_t pthreadAttr;

	// allocate data space
	mPool = new pthread_t[size];
	mThreadData = new threadDataType[size];
	mSlotStack = new int[size];
	mCreateArgs = new createArgumentsType[size];
	mSearchArgs = new searchArgumentsType[size];
	mSearchRangeArgs = new searchRangeArgumentsType[size];

	mPoolSize = size;

	pthread_attr_init(&pthreadAttr);
	pthread_attr_setdetachstate(&pthreadAttr, PTHREAD_CREATE_DETACHED);
	
	for (int i = 0; i < size; i++) {
		// initialize thread data structure
		mThreadData[i].slot = i;
		mThreadData[i].pool = this;
		mThreadData[i].task = 0;
		pthread_mutex_init(&(mThreadData[i].mutex), NULL);
		pthread_cond_init(&(mThreadData[i].condition), NULL);

		// start thread
		int rc = pthread_create(&mPool[i], &pthreadAttr, threadWrapper, (void*) (&(mThreadData[i])));

		rc = rc; // suppress warning in non-debug-mode
		assert(rc == 0);
		
		// put all threads on thread stack
		mSlotStack[i] = i;
	}

	mSlotStackPosition = 0;
	
	pthread_mutex_init(&mSlotStackMutex, NULL);
	pthread_cond_init(&mSlotStackCondition, NULL);
}

// destructor
// stop all threads and discard allocated resources
ThreadPool::~ThreadPool() {
	// wait for locked threads
	wait();

	// lock all threads
	for (int i = 0; i < mPoolSize; i++) {
		getSlot();
	}

	// send task "stop" = 3 to each thread
	for (int i = 0; i < mPoolSize; i++) {
		startSlot(3, i);
	}
	
	// wait for all threads to complete
	wait();

	// release resources
	delete[] mPool;
	delete[] mSlotStack;
	delete[] mCreateArgs;
	delete[] mSearchArgs;
}

// lock next free thread by getting its corresponding
// slot from the slot stack
int ThreadPool::getSlot() {
	int slot;

	// wait for a free slot
	pthread_mutex_lock(&mSlotStackMutex);

	while (mSlotStackPosition >= mPoolSize) {
		pthread_cond_wait(&mSlotStackCondition, &mSlotStackMutex);
	}

	// get next free slot from the stack
	slot = mSlotStack[mSlotStackPosition];
	mSlotStackPosition++;
	
	pthread_mutex_unlock(&mSlotStackMutex);

	return slot;
}

// unlock thread after completing a task by putting
// its corresponding slot back on the slot stack
void ThreadPool::releaseSlot(int slot) {
	pthread_mutex_lock(&mSlotStackMutex);

	mSlotStackPosition--;
	assert(mSlotStackPosition >= 0);
	mSlotStack[mSlotStackPosition] = slot;
	mThreadData[slot].task = 0;	// set task to "wait" = 0

	// signal main thread that "mSlotStackPosition" has changed
	pthread_cond_signal(&mSlotStackCondition);
	
	pthread_mutex_unlock(&mSlotStackMutex);
}

// start locked thread to perform task
void ThreadPool::startSlot(int task, int slot) {
	pthread_mutex_lock(&(mThreadData[slot].mutex));
	
	// set task
	mThreadData[slot].task = task;
	
	// signal thread that "task" has changed
	pthread_cond_signal(&(mThreadData[slot].condition));
	pthread_mutex_unlock(&(mThreadData[slot].mutex));
}

// dispatch a task to create a new MultibitTree
void ThreadPool::createMultibitTree(MultibitTree **tree, Fingerprint **prints, int leafStart, int leafEnd, int nBits, int cardinality, int leafLimit) {
	int slot;

	// lock thread
	slot = getSlot();

	// set attributes
	mCreateArgs[slot].tree = tree;
	mCreateArgs[slot].prints = prints;
	mCreateArgs[slot].leafStart = leafStart;
	mCreateArgs[slot].leafEnd = leafEnd;
	mCreateArgs[slot].nBits = nBits;
	mCreateArgs[slot].cardinality = cardinality;
	mCreateArgs[slot].leafLimit = leafLimit;

	// start thread with task "create" = 1
	startSlot(1, slot);
}

// dispatch a task to search in a MultibitTree
void ThreadPool::searchMultibitTree(MultibitTree *tree, QueryResult *result, Fingerprint *query, int cardinality, float minTanimoto) {
	int slot;

	// lock thread
	slot = getSlot();

	// set attributes	
	mSearchArgs[slot].tree = tree;
	mSearchArgs[slot].result = result;
	mSearchArgs[slot].query = query;
	mSearchArgs[slot].cardinality = cardinality;
	mSearchArgs[slot].minTanimoto = minTanimoto;

	// start thread with task "search" = 2
	startSlot(2, slot);
}

// dispatch a task to search in an array of MultibitTrees
void ThreadPool::searchMultibitTreeRange(MultibitTree **trees, int min, int max, QueryResult *result, Fingerprint *query, int cardinality, float minTanimoto) {
	int slot;

	// lock thread
	slot = getSlot();

	// set attributes	
	mSearchRangeArgs[slot].trees = trees;
	mSearchRangeArgs[slot].min = min;
	mSearchRangeArgs[slot].max = max;
	mSearchRangeArgs[slot].result = result;
	mSearchRangeArgs[slot].query = query;
	mSearchRangeArgs[slot].cardinality = cardinality;
	mSearchRangeArgs[slot].minTanimoto = minTanimoto;

	// start thread with task "searchRange" = 4
	startSlot(4, slot);
}

// wait until all threads have completed
void ThreadPool::wait() {
	pthread_mutex_lock(&mSlotStackMutex);

	while (mSlotStackPosition > 0) {
		pthread_cond_wait(&mSlotStackCondition, &mSlotStackMutex);
	}

	pthread_mutex_unlock(&mSlotStackMutex);
}


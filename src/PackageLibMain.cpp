// PackageLibMain.cpp
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

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <string.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "Misc.h"
#include "Grid1D.h"

// This file contains the pure c-functions for the R-library-interface.
// R_init_useCall	register .Call-Methods
// mbtLoadCall		wrapper for Grid1D-constructor
// mbtSearchCall	wrapper for Grid1D::search
// mbtSearchFileCall	wrapper for Grid1D::searchFile
// mbtUnloadCall	wrapper for Grid1D-destructor
// mbtStatistics	wrapper for 

// maximal line size = length of ascii representation of fingerprint
#define STRSIZE 4000

#ifdef __cplusplus
extern "C" {
#endif

static Grid1D *grid = NULL;		// static Grid1D as data structure singleton

// call destructor of static grid datastructure
void mbtUnload() {
	if (grid != NULL) {
		delete grid;
		grid = NULL;
	}
}

// check, if a character is considered to be a white-space or seperator
int isWS(char c) {
	return (
		(c == '"') ||
		(c == '\'') ||
		(c == ',') ||
		(c == ';') ||
		(c == ' ') ||
		(c == '\t')
	);
}

// check, if a character is considered to be end of line
int isEOL(char c) {
	return (
		(c == 10) ||
		(c == 13) ||
		(c == 0)
	);
}

// parse prints from file, return number of parsed fields
int parseLine(FILE *in, char *str, long long *idx1, long long *end1, long long *idx2, long long *end2) {
	// read line from file
	if (fgets(str, STRSIZE, in) == NULL) {
		return 0;
	}
	
	// parse line
	
	// find start of first field
	*idx1 = 0;
	while (isWS(str[*idx1]) && (*idx1 < STRSIZE-1)) {
		(*idx1)++;
	}

	// find end of first field
	*end1 = *idx1;
	while (!isWS(str[*end1]) && (*end1 < STRSIZE-1) && !isEOL(str[*end1])) {
		(*end1)++;
	}

	if (!isEOL(str[*end1])) {
		// find start of second field
		*idx2 = *end1;
		while (isWS(str[*idx2]) && (*idx2 < STRSIZE-1)) {
			(*idx2)++;
		}
		
		// find end of second field
		*end2 = *idx2;
		while (!isWS(str[*end2]) && (*end2 < STRSIZE-1) && !isEOL(str[*end2])) {
			(*end2)++;
		}
	} else {
		*idx2 = 0;
		*end2 = 0;
	}

	// cut field 1
	str[*end1] = 0;

	if (*idx2 != *end2) {
		// if there are two fields, cut field 2
		str[*end2] = 0;
		return 2;
	}

	return 1;
}

// read input file and call constructor for static data structure
int mbtLoad(const char *filename, int threads, long long size, int leafLimit) {
	Fingerprint **prints;
	long long sizePrints;
	int nBits;
	int len;
	char str[STRSIZE];
	long long fields;
	long long idx1, end1, idx2, end2;
	char *idStr;
	FILE *in;

	// initialize Fingerprint data structure (cardinality-map)
        Fingerprint::init();

	// delete an existing grid
	if (grid != NULL) {
		mbtUnload();
	}

	sizePrints = size;

	// count number of prints in file
	if (sizePrints == 0) {
		in = fopen(filename, "r");

		if (in == NULL) {
			return(0);
		}

                while (fgets(str, STRSIZE, in)) {
			sizePrints++;
		}

		fclose(in);
	}

	nBits = 0;
        in = fopen(filename, "r");

	if (in == NULL) {
		return(0);
	}

	// create array of Fingerprints
        prints = new Fingerprint*[sizePrints];

	// read prints from file
        for (long long i = 0; i < sizePrints; i++) {
		// parse line from file
		fields = parseLine(in, str, &idx1, &end1, &idx2, &end2);

		if (fields == 0) {
			sizePrints = i;
			break;
                }
		
		if (fields == 1) {
			// if there is only one field, use line as id
			idStr = new char[13];
			sprintf(idStr, "%012" PRId64, i+1);
			// create fingerprint
			prints[i] = new Fingerprint(idStr, str+idx1);
		} else {
			// if there are two fields, use first string as id
			idStr = new char[end1-idx1+1];
			strcpy(idStr, str+idx1);
			// create fingerprint
			prints[i] = new Fingerprint(idStr, str+idx2);
		}

		// compute maximal length of Fingerprints
                len = prints[i]->getLength();
                if (len > nBits) {
                        nBits = len;
                }
        }

        fclose(in);

	// build Grid1D data structure
	grid = new Grid1D(prints, sizePrints, nBits, threads, leafLimit);

	return(sizePrints);
}

// recursively copy QueryResults into R vectors for prints and tanimoto coefficients
void insertQueryResultNodesRec(long long *idx, SEXP prints, double *tanimotosPtr, QueryResultNode *node, long long sizeResult) {
	if ((node != NULL) && ((*idx) < sizeResult)) {
		insertQueryResultNodesRec(idx, prints, tanimotosPtr, node->mLeft, sizeResult);

		SET_STRING_ELT(prints, *idx, mkChar(node->mPrint->getId()));
		tanimotosPtr[*idx] = node->mTanimoto;		// copy tanimoto coefficent

		(*idx)++;

		insertQueryResultNodesRec(idx, prints, tanimotosPtr, node->mRight, sizeResult);
	}
}

// copy QueryResults into R vectors for prints and tanimoto coefficients
void insertQueryResultNodes(long long *idx, SEXP prints, double *tanimotosPtr, QueryResultNode *node, long long sizeResult) {
	while ((node != NULL) && ((*idx) < sizeResult)) {

		SET_STRING_ELT(prints, *idx, mkChar(node->mPrint->getId()));
		tanimotosPtr[*idx] = node->mTanimoto;		// copy tanimoto coefficent

		(*idx)++;
		node = node->mRight;
	}
}

// copy QueryResults into R vectors for queries, prints and tanimoto coefficients
void insertQueryResultNodesWithId(long long *idx, SEXP queries, SEXP prints, double *tanimotosPtr, QueryResultNode *node, long long sizeResult) {
	while ((node != NULL) && ((*idx) < sizeResult)) {
		SET_STRING_ELT(queries, *idx, mkChar(node->mQueryId));
		SET_STRING_ELT(prints, *idx, mkChar(node->mPrint->getId()));
		tanimotosPtr[*idx] = node->mTanimoto;		// copy tanimoto coefficent

		(*idx)++;
		node = node->mRight;
	}
}

// call Grid1D::search and store results into vector of vectors
SEXP mbtSearch(const char *query, double minTanimoto, long long size, int sort) {
	SEXP result;
	SEXP names;
	SEXP prints;
	SEXP tanimotos;
	double *tanimotosPtr;
	QueryResult queryResult(sort, NULL, NULL);
	Fingerprint queryPrint(NULL, query);
	long long idx;
	long long sizeResult;

	// call search-method
	if (grid != NULL) {
		grid->initStatistics();
		grid->search(&queryResult, &queryPrint, minTanimoto);
		grid->setSizeLastSearch(1);
	}

	sizeResult = queryResult.getSize();

	if (size > 0) {
		sizeResult = MIN(size, sizeResult);
	}

	// allocate R data structures for result
	PROTECT(prints = allocVector(STRSXP, sizeResult));
	PROTECT(tanimotos = allocVector(REALSXP, sizeResult));
	tanimotosPtr = REAL(tanimotos);

	// copy result into R data structures
	idx = 0;
	if (sort) {
		insertQueryResultNodesRec(&idx, prints, tanimotosPtr, queryResult.getRootNode(), sizeResult);
	} else {
		insertQueryResultNodes(&idx, prints, tanimotosPtr, queryResult.getRootNode(), sizeResult);
	}

	// allocate vector for the two result vectors
	PROTECT(result = allocVector(VECSXP, 2));

	SET_VECTOR_ELT(result, 0, prints);
	SET_VECTOR_ELT(result, 1, tanimotos);

	// set name attributes for the two result vectors
	PROTECT(names = allocVector(STRSXP, 2));

	SET_STRING_ELT(names, 0, mkChar("fingerprint"));
	SET_STRING_ELT(names, 1, mkChar("tanimoto"));
	setAttrib(result, R_NamesSymbol, names);

	UNPROTECT(4);

	return(result);
}

// call Grid1D::search for each fingerprint in file and store results into vector of vectors
// if a result file is specified, write the results in to this file and return nothing to the R-function
SEXP mbtSearchFile(const char *filename, double minTanimoto, const char *resultFile, const char *seperator) {
	SEXP result;
	SEXP names;
	SEXP queries;
	SEXP prints;
	SEXP tanimotos;
	double *tanimotosPtr;
	Fingerprint *queryPrint;
	long long idx;
	long long sizeResult;
	FILE *in;
	FILE *out = NULL;
	char str[STRSIZE];
	long long fields;
	long long idx1, end1, idx2, end2;
	long long i;
	char *idStr;

	if ((resultFile != NULL) && (resultFile[0] != 0) && (seperator != NULL)) {
		// if specified, open result file
		out = fopen(resultFile, "w");
		// print column headers
		fprintf(out, "query%sfingerprint%stanimoto\n", seperator, seperator);
	}

	QueryResult queryResult(0, out, seperator);

	if (grid != NULL) {
		grid->initStatistics();

		// open input file
		in = fopen(filename, "r");

		if (in != NULL) {
			i = 0;
			while (1) {
				// for each line parse fingerprint
				fields = parseLine(in, str, &idx1, &end1, &idx2, &end2);
				if (fields == 0) {
					break;
				}
				
				// create query fingerprint
				if (fields == 1) {
					// if there is only one field, use line as id
					idStr = new char[13];
					sprintf(idStr, "%012" PRId64, i+1);
					// create fingerprint
					queryPrint = new Fingerprint(idStr, str+idx1);
				} else {
					// if there are two fields, use first string as id
					idStr = new char[end1-idx1+1];
					strcpy(idStr, str+idx1);
					// create fingerprint
					queryPrint = new Fingerprint(idStr, str+idx2);
				}

				// call asychonous search-method
				grid->searchAsync(&queryResult, queryPrint, minTanimoto);
				i++;
			}

			grid->setSizeLastSearch(i);
      
			// wait for running threads
			grid->wait();
			fclose(in);
		}
	}

	sizeResult = queryResult.getSize();

	if (out != NULL) {
		// if a result file was specified close it and return NULL
		fclose(out);
		return(R_NilValue);
	}

	// allocate R data structures for result
	PROTECT(queries = allocVector(STRSXP, sizeResult));
	PROTECT(prints = allocVector(STRSXP, sizeResult));
	PROTECT(tanimotos = allocVector(REALSXP, sizeResult));
	tanimotosPtr = REAL(tanimotos);

	// copy result into R data structures
	idx = 0;
	insertQueryResultNodesWithId(&idx, queries, prints, tanimotosPtr, queryResult.getRootNode(), sizeResult);

	// allocate vector for the three result vectors
	PROTECT(result = allocVector(VECSXP, 3));

	SET_VECTOR_ELT(result, 0, queries);
	SET_VECTOR_ELT(result, 1, prints);
	SET_VECTOR_ELT(result, 2, tanimotos);

	// set name attributes for the two result vectors
	PROTECT(names = allocVector(STRSXP, 3));

	SET_STRING_ELT(names, 0, mkChar("query"));
	SET_STRING_ELT(names, 1, mkChar("fingerprint"));
	SET_STRING_ELT(names, 2, mkChar("tanimoto"));
	setAttrib(result, R_NamesSymbol, names);

	UNPROTECT(5);

	return(result);
}

// call Grid1D::getStatistics
SEXP mbtStatistics() {
	SEXP result;
	SEXP names;
	SEXP values;
	SEXP params;
	SEXP percents;
	
	double *valuesPtr;
	double *percentsPtr;

	PROTECT(params = allocVector(STRSXP, 3));
	PROTECT(values = allocVector(REALSXP, 3));
	PROTECT(percents = allocVector(REALSXP, 3));

	valuesPtr = REAL(values);
	percentsPtr = REAL(percents);

	if (grid != NULL) {
		grid->getStatistics(valuesPtr, percentsPtr);
	}

	SET_STRING_ELT(params, 0, mkChar("XOR-Hash"));
	SET_STRING_ELT(params, 1, mkChar("Tanimoto"));
	SET_STRING_ELT(params, 2, mkChar("Total"));

	PROTECT(result = allocVector(VECSXP, 3));

	SET_VECTOR_ELT(result, 0, params);
	SET_VECTOR_ELT(result, 1, values);
	SET_VECTOR_ELT(result, 2, percents);
	
	// set name attributes for the two result vectors
	PROTECT(names = allocVector(STRSXP, 3));

	SET_STRING_ELT(names, 0, mkChar("Checkpoint"));
	SET_STRING_ELT(names, 1, mkChar("Count"));
	SET_STRING_ELT(names, 2, mkChar("Percentage"));
	setAttrib(result, R_NamesSymbol, names);

 	UNPROTECT(5);

	return(result);
}

// wrapper for R-function mbtLoadCall
SEXP mbtLoadCall(SEXP filename, SEXP threads, SEXP size, SEXP leafLimit) {
	SEXP result;
	int *resultPtr;

	PROTECT(filename = AS_CHARACTER(filename));
	PROTECT(threads = AS_INTEGER(threads));
	PROTECT(size = AS_INTEGER(size));
	PROTECT(leafLimit = AS_INTEGER(leafLimit));

	PROTECT(result = NEW_INTEGER(1));
	resultPtr = INTEGER_POINTER(result);

	resultPtr[0] = mbtLoad(CHAR(STRING_ELT(filename, 0)), INTEGER_POINTER(threads)[0], INTEGER_POINTER(size)[0], INTEGER_POINTER(leafLimit)[0]);

	UNPROTECT(5);

	return(result);
}

// wrapper for R-function mbtSearchCall
SEXP mbtSearchCall(SEXP query, SEXP minTanimoto, SEXP size, SEXP sort) {
	SEXP result;

	PROTECT(query = AS_CHARACTER(query));
	PROTECT(minTanimoto = AS_NUMERIC(minTanimoto));
	PROTECT(size = AS_INTEGER(size));
	PROTECT(sort = AS_INTEGER(sort));

	result = mbtSearch(CHAR(STRING_ELT(query, 0)), REAL(minTanimoto)[0], INTEGER_POINTER(size)[0], INTEGER_POINTER(sort)[0]);
	
	UNPROTECT(4);

	return(result);
}

// wrapper for R-function mbtSearchFileCall
SEXP mbtSearchFileCall(SEXP filename, SEXP minTanimoto, SEXP resultFile, SEXP seperator) {
	SEXP result;

	PROTECT(filename = AS_CHARACTER(filename));
	PROTECT(minTanimoto = AS_NUMERIC(minTanimoto));
	PROTECT(resultFile = AS_CHARACTER(resultFile));
	PROTECT(seperator = AS_CHARACTER(seperator));

	result = mbtSearchFile(CHAR(STRING_ELT(filename, 0)), REAL(minTanimoto)[0], CHAR(STRING_ELT(resultFile, 0)), CHAR(STRING_ELT(seperator, 0)));
	
	UNPROTECT(4);

	return(result);
}

// wrapper for R-function mbtUnloadCall
SEXP mbtUnloadCall() {
	mbtUnload();

	return(R_NilValue);
}

SEXP mbtStatisticsCall() {
	SEXP result;

	result = mbtStatistics();
  
	return(result);
}

// register wrapper-functions
void R_init_useCall(DllInfo *info) {
	R_CallMethodDef callMethods[]  = {
	  {"mbtLoadCall", (DL_FUNC) &mbtLoadCall, 4},
	  {"mbtSearchCall", (DL_FUNC) &mbtSearchCall, 4},
	  {"mbtSearchFileCall", (DL_FUNC) &mbtSearchFileCall, 4},
	  {"mbtUnloadCall", (DL_FUNC) &mbtUnloadCall, 0},
	  {"mbtStatisticsCall", (DL_FUNC) &mbtStatisticsCall, 0},
	  {NULL, NULL, 0}
	};
	
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

#ifdef __cplusplus
}
#endif

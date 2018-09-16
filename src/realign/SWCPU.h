/*
 * SWCPU.h
 *
 *  Created on: Jun 15, 2011
 *      Author: fritz
 */

#ifndef SWCPU_H_
#define SWCPU_H_

#define pRef pBuffer1
#define pQry pBuffer2

#include "IAlignment.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using std::endl;
using std::cout;
using std::max;

#define CIGAR_STOP 10
#define short_min -16000
#define result_number 4
#define line_end '\0'
#define ref_position 0
#define qstart 1
#define qend 2
#define alignment_offset 3
#define param_best_read_index 0
#define param_best_ref_index 1

#define CIGAR_M 0
#define CIGAR_I 1
#define CIGAR_D 2
#define CIGAR_N 3
#define CIGAR_S 4
#define CIGAR_H 5
#define CIGAR_P 6
#define CIGAR_EQ 7
#define CIGAR_X 8


typedef float Score;

struct MatrixElement {
	Score score;
	int indelRun;
	char direction;
};


const char trans[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0,
		4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4,
		4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

class SWCPUCor: public IAlignment {
public:
	SWCPUCor(int gpu_id);
	virtual ~SWCPUCor();
	virtual int GetScoreBatchSize() const;
	virtual int GetAlignBatchSize() const;
	virtual int BatchScore(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, float * const results,
			void * extData);
	virtual int BatchAlign(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, Align * const results,
			void * extData);

	virtual int SingleAlign(int const mode, int const corridor,
			char const * const refSeq, char const * const qrySeq,
			Align & result, void * extData);
private:

	//bool cigar;
	//short scores[6][6];
	Score mat;
	Score mis;
	Score gap_open_read;
	Score gap_open_ref;
	Score gap_ext;
	Score gap_ext_min;
	Score gap_decay;

	MatrixElement * alignMatrix;
	int * binaryCigar;

	//meta info
	unsigned int batch_size; //effictive thread number that is started per call

	int printCigarElement(char const op, int const length, char * cigar);

	int computeCigarMD(Align & result, int const gpuCigarOffset,
			int const * const gpuCigar, char const * const refSeq, int corr_length, int read_length, int const QStart, int const QEnd);

	Score SW_Score(char const * const scaff, char const * const read, int * result, int corr_length, MatrixElement * mat_pointer);

	bool Backtracking_CIGAR(char const * const scaff, char const * const read,
			int *& result, int *& alignments, int corr_length, int read_length, int alignment_length, MatrixElement * mat_pointer);

	void print_matrix(int alignment_length, const char* const refSeq,
			int read_length, const char* const qrySeq, int corr_length,
			MatrixElement* mat_pointer);
};

#endif /* SWCPU_H_ */

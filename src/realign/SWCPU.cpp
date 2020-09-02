/*
 * SWCPU.cpp
 *
 *  Created on: Jun 15, 2011
 *      Author: fritz
 */

#include "SWCPU.h"

//TODO: hack to pass data for debug output
Align cur_align;

SWCPUCor::SWCPUCor(int gpu_id) {

//	cigar = bool(((gpu_id >> 8) & 0xFF) == 1);

	batch_size = 1;

	mat = 2.0f;
	mis = -5.0f;
	gap_open_read = -5.0f;
	gap_open_ref = -5.0f;
	gap_ext = -5.0f;
	gap_decay = 0.05f;
	gap_ext_min = -1.0f;

	long maxLen = (long) 100000 * (long) 20000;
	alignMatrix = new MatrixElement[maxLen];

	fprintf(stderr, "Allocationg: %lu\n", maxLen * sizeof(MatrixElement));

	binaryCigar = new int[200000];

//	short temp[6][6] = { mat, mis, mis, mis, 0, mis, mis, mat, mis, mis, 0, mis,
//			mis, mis, mat, mis, 0, mis, mis, mis, mis, mat, 0, mis, 0, 0, 0, 0,
//			0, 0, 0, 0, 0, 0, 0, mat };
//	memcpy(scores, temp, 6 * 6 * sizeof(short));

	fprintf(stderr, "SWCPU initialized\n");
}

SWCPUCor::~SWCPUCor() {
	delete[] alignMatrix;
	alignMatrix = 0;
	delete[] binaryCigar;
	binaryCigar = 0;
}

Score SWCPUCor::SW_Score(char const * const refSeqList,
		char const * const qrySeqList, int * fwResults, int corr_length,
		MatrixElement * mat_pointer) {

//	memset(local_mat_line, 0, corr_length * sizeof(short));
	char const * scaff = refSeqList;
	char const * read = qrySeqList;

	//Init matrix lines
	MatrixElement * matrix = mat_pointer;

	//Init matrix lines
	for (int i = 0; i < corr_length; ++i) {
		//local_mat_line[i] = 0;
		matrix[i].direction = CIGAR_STOP;
		matrix[i].indelRun = 0;
		matrix[i].score = 0;
	}
	matrix[corr_length].direction = CIGAR_STOP;
	matrix[corr_length].indelRun = 0;
	matrix[corr_length].score = 0;

	Score curr_max = -1.0f;
	int read_index = 0;

	int x = 0;
	for (; *read != line_end; ++read) {
		char read_char_cache = *read;
		matrix += (corr_length + 1);
//		short left_cell = 0;
		matrix[0].direction = CIGAR_STOP;
		matrix[0].indelRun = 0;
		matrix[0].score = 0;

		for (int ref_index = 0; ref_index < corr_length - 1; ++ref_index) {

			MatrixElement & diag = matrix[-(corr_length + 1) + ref_index + 1];
			MatrixElement & up = matrix[-(corr_length + 1) + ref_index + 2];
			MatrixElement & left = matrix[ref_index];

			bool eq = read_char_cache == scaff[ref_index];
			Score diag_cell = diag.score + ((eq) ? mat : mis);

			Score up_cell = 0;
			Score left_cell = 0;

			int ins_run = 0;
			int del_run = 0;

			if (up.direction == CIGAR_I) {
				ins_run = up.indelRun;
				if (up.score == 0) {
					up_cell = 0;
				} else {
					up_cell = up.score
							+ std::min(gap_ext_min,
									gap_ext + ins_run * gap_decay);
				}
			} else {
				up_cell = up.score + gap_open_read;
			}

			if (left.direction == CIGAR_D) {
				del_run = left.indelRun;
				if (left.score == 0) {
					left_cell = 0;
				} else {
					left_cell = left.score
							+ std::min(gap_ext_min,
									gap_ext + del_run * gap_decay);
				}
			} else {
				left_cell = left.score + gap_open_ref;
			}

			//find max
			Score max_cell = 0;
			max_cell = max(left_cell, max_cell);
			max_cell = max(diag_cell, max_cell);
			max_cell = max(up_cell, max_cell);

			MatrixElement & current = matrix[(ref_index + 1)];
			if (del_run > 0 && max_cell == left_cell) {
				current.score = max_cell;
				current.direction = CIGAR_D;
				current.indelRun = del_run + 1;
			} else if (ins_run > 0 && max_cell == up_cell) {
				current.score = max_cell;
				current.direction = CIGAR_I;
				current.indelRun = ins_run + 1;
			} else if (max_cell == diag_cell) {
				current.score = max_cell;
				if (eq) {
					current.direction = CIGAR_EQ;
				} else {
					current.direction = CIGAR_X;
				}
				current.indelRun = 0;
			} else if (max_cell == left_cell) {
				current.score = max_cell;
				current.direction = CIGAR_D;
				current.indelRun = 1;
			} else if (max_cell == up_cell) {
				current.score = max_cell;
				current.direction = CIGAR_I;
				current.indelRun = 1;
			} else {
				current.score = 0;
				current.direction = CIGAR_STOP;
				current.indelRun = 0;
			}

			if (max_cell > curr_max) {
				curr_max = max_cell;
				fwResults[param_best_ref_index] = ref_index;
				fwResults[param_best_read_index] = read_index;
				fwResults[3] = curr_max;

			}
		}
		matrix[corr_length].direction = CIGAR_STOP;
		matrix[corr_length].score = 0;
		matrix[corr_length].indelRun = 0;

		scaff++;
		read_index += 1;
	}
	fwResults[2] = (read_index - fwResults[0]) - 1;
	if (read_index == 0) {
		fwResults[0] = fwResults[1] = 2;
	}
	return curr_max;
}

int SWCPUCor::printCigarElement(char const op, int const length, char * cigar) {
	int offset = 0;
	offset = sprintf(cigar, "%d%c", length, op);
	return offset;
}

int SWCPUCor::computeCigarMD(Align & result, int const gpuCigarOffset,
		int const * const gpuCigar, char const * const refSeq, int corr_length,
		int read_length, int const QStart, int const QEnd) {
	int alignment_length = corr_length + read_length + 1;

	int finalCigarLength = 0;

	int cigar_offset = 0;
//	int md_offset = 0;

	if (((gpuCigar[gpuCigarOffset] >> 4) + QStart) > 0) {
		fprintf(stderr, "Adding %d to QSTart\n", QStart);
		result.QStart = (gpuCigar[gpuCigarOffset] >> 4) + QStart;
		cigar_offset += printCigarElement('S', result.QStart,
				result.pRef + cigar_offset);
		finalCigarLength += result.QStart;
	}

	int cigar_m_length = 0;
//	int md_eq_length = 0;
	int ref_index = 0;
	for (int j = gpuCigarOffset + 1; j < (alignment_length - 1); ++j) {
		int op = gpuCigar[j] & 15;
		int length = gpuCigar[j] >> 4;

		//debugCigar(op, length);

		switch (op) {
		case CIGAR_X:
			cigar_m_length += length;

			//Produces: 	[0-9]+(([A-Z]+|\^[A-Z]+)[0-9]+)*
			//instead of: 	[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
//			md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
//			for (int k = 0; k < length; ++k) {
//				md_offset += sprintf(result.pQry + md_offset, "%c",
//						refSeq[ref_index++]);
//			}
//			md_eq_length = 0;

			break;
		case CIGAR_EQ:
			cigar_m_length += length;
//			md_eq_length += length;
			ref_index += length;
			break;
		case CIGAR_D:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length,
						result.pRef + cigar_offset);
				finalCigarLength += cigar_m_length;
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('D', length,
					result.pRef + cigar_offset);

//			md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
//			md_eq_length = 0;
//			result.pQry[md_offset++] = '^';
//			for (int k = 0; k < length; ++k) {
//				result.pQry[md_offset++] = refSeq[ref_index++];
//			}

			break;
		case CIGAR_I:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length,
						result.pRef + cigar_offset);
				finalCigarLength += cigar_m_length;
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('I', length,
					result.pRef + cigar_offset);
			finalCigarLength += length;

			break;
		default:
			fprintf(stderr, "Invalid cigar string: %d\n", op);
			std::cout << "Offset: " << gpuCigarOffset << std::endl;
			for (int x = 0; x < alignment_length * 2; ++x) {
				std::cout << gpuCigar[x] << " ";
			}
			std::cout << std::endl;
			exit(1);
		}
	}
//	md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
	if (cigar_m_length > 0) {
		cigar_offset += printCigarElement('M', cigar_m_length,
				result.pRef + cigar_offset);
		finalCigarLength += cigar_m_length;
		cigar_m_length = 0;
	}

	if (((gpuCigar[alignment_length - 1] >> 4) + QEnd) > 0) {
		fprintf(stderr, "Adding %d to QEnd\n", QEnd);
		result.QEnd = (gpuCigar[alignment_length - 1] >> 4) + QEnd;
		cigar_offset += printCigarElement('S', result.QEnd,
				result.pRef + cigar_offset);
		finalCigarLength += result.QEnd;

	}
	//TODO: fix
	result.Identity = 1.0f;
	result.pRef[cigar_offset] = '\0';
//	result.pQry[md_offset] = '\0';

	return finalCigarLength;
}

bool SWCPUCor::Backtracking_CIGAR(char const * const scaff,
		char const * const read, int *& fwdResults, int *& alignments,
		int corr_length, int read_length, int alignment_length,
		MatrixElement * mat_pointer) {

	bool valid = true;

	MatrixElement * matrix = mat_pointer;

	int best_read_index = fwdResults[param_best_read_index];
	int best_ref_index = fwdResults[param_best_ref_index];

	int cigarLenth = 0;
	int cigarLengthCheck = 0;

	int totalDelLength = 0;
	int totalINsLength = 0;

	int minCorridor = corr_length * 0.01f;
	int maxCorridor = corr_length - minCorridor;

	if (best_read_index > 0) {
		matrix += (((corr_length + 1) * (best_read_index + 1)));

		int abs_ref_index = best_ref_index + best_read_index;
		int alignment_index = alignment_length - 1;

		int pointer = CIGAR_STOP;
		int cigar_element = CIGAR_S;
		int cigar_length = fwdResults[qend];
		cigarLenth += fwdResults[qend];
		while ((pointer = matrix[(best_ref_index + 1)].direction) != CIGAR_STOP) {
//			Log.Message("Best ref index: %d (%d)", best_ref_index + 1, (corr_length + 1));
//			printf("%s\t%d\t%d\t%d\t%d\n", (char *) cur_align.ExtendedData,
//					cur_align.NM, best_read_index, best_ref_index + 1,
//					corr_length + 1);
			if (best_ref_index <= minCorridor
					|| best_ref_index >= maxCorridor) {
				fprintf(stderr, "Corridor probably too small\n");
				valid = false;
//				getchar();
			}

			if (pointer == CIGAR_X || pointer == CIGAR_EQ) {
				matrix -= ((corr_length + 1));
				best_read_index -= 1;
				abs_ref_index -= 1;

				cigarLenth += 1;
			} else if (pointer == CIGAR_I) {
				matrix -= ((corr_length + 1));
				best_read_index -= 1;
				best_ref_index += 1;

				cigarLenth += 1;
			} else if (pointer == CIGAR_D) {
				best_ref_index -= 1;
				abs_ref_index -= 1;

			} else {
				fprintf(stderr, "Error in backtracking. Invalid CIGAR operation found\n");
				exit(1);

			}

			if (pointer == cigar_element) {
				cigar_length += 1;
			} else {
				alignments[alignment_index--] = (cigar_length << 4
						| cigar_element);
				if (cigar_element != CIGAR_D) {
					cigarLengthCheck += cigar_length;
				}

				cigar_element = pointer;
				cigar_length = 1;
			}
		}
		alignments[alignment_index--] = (cigar_length << 4 | cigar_element);
		if (cigar_element != CIGAR_D) {
			cigarLengthCheck += cigar_length;
		}

		alignments[alignment_index] = ((best_read_index + 1) << 4 | CIGAR_S);
		cigarLengthCheck += (best_read_index + 1);
		cigarLenth += (best_read_index + 1);
		fwdResults[ref_position] = abs_ref_index + 1;
		fwdResults[qstart] = best_read_index + 1;
		//qend was set by "forward" kernel
		fwdResults[alignment_offset] = alignment_index;

		if (cigarLenth != cigarLengthCheck) {
			fprintf(stderr, "Error in CIGAR length: %d vs %d\n", cigarLenth, cigarLengthCheck);
		} else {
			if (read_length != cigarLenth) {
				fprintf(stderr, "Error read length != cigar length: %d vs %d\n", read_length, cigarLenth);
				exit(1);
			}
		}
		fprintf(stderr, "Read length: %d, CIGAR length: %d\n", read_length, cigarLenth);
	}
	return valid;
}

int SWCPUCor::GetScoreBatchSize() const {
	return 0;
}
int SWCPUCor::GetAlignBatchSize() const {
	return 0;
}

int SWCPUCor::BatchAlign(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, Align * const results,
		void * extData) {

	throw "Not implemented";

	fprintf(stderr, "Unsupported alignment mode %i\n", mode);
	return 0;
}

void SWCPUCor::print_matrix(int alignment_length, const char* const refSeq,
		int read_length, const char* const qrySeq, int corr_length,
		MatrixElement* mat_pointer) {

	printf("    - ");
	for (int x = 0; x < alignment_length - 1; ++x) {
		printf(" %c ", refSeq[x]);
	}
	printf("\n");
	for (size_t row = 0; row < read_length + 1; ++row) {
		if (row == 0) {
			printf("-: ");
		} else {
			printf("%c: ", qrySeq[row - 1]);
		}
		for (int x = 0; x < row; ++x) {
			printf("   ");
		}
		for (size_t col = 0; col < corr_length + 1; ++col) {
			MatrixElement* cell = mat_pointer + (row * (corr_length + 1) + col);
			printf("%*d ", 2, cell->indelRun);
		}
		printf("\n");
	}

	printf("    - ");
	for (int x = 0; x < alignment_length - 1; ++x) {
		printf(" %c ", refSeq[x]);
	}
	printf("\n");
	for (size_t row = 0; row < read_length + 1; ++row) {
		if (row == 0) {
			printf("-: ");
		} else {
			printf("%c: ", qrySeq[row - 1]);
		}
		for (int x = 0; x < row; ++x) {
			printf("   ");
		}
		for (size_t col = 0; col < corr_length + 1; ++col) {
			MatrixElement* cell = mat_pointer + (row * (corr_length + 1) + col);
			printf("%*d ", 2, cell->direction);
		}
		printf("\n");
	}

	printf("    - ");
	for (int x = 0; x < alignment_length - 1; ++x) {
		printf(" %c ", refSeq[x]);
	}
	printf("\n");
	for (size_t row = 0; row < read_length + 1; ++row) {
		if (row == 0) {
			printf("-: ");
		} else {
			printf("%c: ", qrySeq[row - 1]);
		}
		for (int x = 0; x < row; ++x) {
			printf("   ");
		}
		for (size_t col = 0; col < corr_length + 1; ++col) {
			MatrixElement* cell = mat_pointer + (row * (corr_length + 1) + col);
			printf("%*.*f ", 2, 0, cell->score);
		}
		printf("\n");
	}
}

int SWCPUCor::SingleAlign(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, Align & align,
		void * extData) {

//	Log.Message("Aligning: ");
//	Log.Message("%s", refSeq);
//	Log.Message("%s", qrySeq);

	cur_align = align;

	int * clipping = 0;
	if (extData == 0) {
		clipping = new int[2];
		clipping[0] = 0;
		clipping[1] = 0;
	} else {
		clipping = (int *) extData;
	}

	int read_length = strlen(qrySeq);
	fprintf(stderr, "Read length (single align) is %d\n", read_length);
	align.pBuffer1 = new char[read_length * 4];
	//	align.pBuffer2 = new char[read_length * 4];
	align.pBuffer2 = new char[1];
	align.pBuffer2[0] = '\0';

	int finalCigarLength = 0;

	int corr_length = corridor;
	int alignment_length = (corr_length + read_length + 1);

	int * fwdResults = new int[result_number];

	Score score = SW_Score(refSeq, qrySeq, fwdResults, corr_length,
			alignMatrix);

//	print_matrix(alignment_length, refSeq, read_length, qrySeq, corr_length,
//			alignMatrix);
//	Log.Message("%d, %d, %d, %d", fwdResults[0], fwdResults[1], fwdResults[2], fwdResults[3]);

	bool valid = Backtracking_CIGAR(refSeq, qrySeq, fwdResults, binaryCigar,
			corr_length, read_length, alignment_length, alignMatrix);

	if (valid) {
		finalCigarLength = computeCigarMD(align, fwdResults[3], binaryCigar,
				refSeq + fwdResults[0], corr_length, read_length, clipping[0],
				clipping[1]);
		align.PositionOffset = fwdResults[0];
		align.Score = score;
	}

	delete[] fwdResults;

	if (extData == 0) {
		delete[] clipping;
		clipping = 0;
	}

	if (!valid) {
		finalCigarLength = -1;
	}

	return finalCigarLength;

}

int SWCPUCor::BatchScore(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, float * const results,
		void * extData) {
	throw "Not implemented";
	return 0;

}


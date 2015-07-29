/*
 * SWCPU.cpp
 *
 *  Created on: Jun 15, 2011
 *      Author: fritz
 */

#include "SWCPU.h"

SWCPUCor::SWCPUCor(int gpu_id) {

	cigar = bool(((gpu_id >> 8) & 0xFF) == 1);

	batch_size = 1;

	mat = 1;
	mis = -1;
	gap_read = -1;
	gap_ref = -1;

	short temp[6][6] = { mat, mis, mis, mis, 0, mis, mis, mat, mis, mis, 0, mis,
			mis, mis, mat, mis, 0, mis, mis, mis, mis, mat, 0, mis, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, mat };

	memcpy(scores, temp, 6 * 6 * sizeof(short));

	Log.Message("SWCPU initialized");
}

SWCPUCor::~SWCPUCor() {

}

float SWCPUCor::SW_Score(char const * const refSeqList,
		char const * const qrySeqList, short * fwResults, int corr_length,
		MatrixElement * mat_pointer, short * local_mat_line) {

	memset(local_mat_line, 0, corr_length * sizeof(short));
	char const * scaff = refSeqList;
	char const * read = qrySeqList;

	//Init matrix lines
	MatrixElement * matrix = mat_pointer;

	//Init matrix lines
	for (short i = 0; i < corr_length; ++i) {
		local_mat_line[i] = 0;
		matrix[i].direction = CIGAR_STOP;
		matrix[i].indelRun = 0;
		matrix[i].score = 0;
	}
	matrix[corr_length].direction = CIGAR_STOP;
	matrix[corr_length].indelRun = 0;
	matrix[corr_length].score = 0;

	short curr_max = -1;
	short read_index = 0;

	int x = 0;
	for (; *read != line_end; ++read) {
		char read_char_cache = *read;
		matrix += (corr_length + 1);
		short left_cell = 0;
		matrix[0].direction = CIGAR_STOP;
		matrix[0].indelRun = 0;
		matrix[0].score = 0;

		for (short ref_index = 0; ref_index < corr_length - 1; ++ref_index) {

			//init values
			left_cell += gap_ref;
			short diag_cell = local_mat_line[ref_index];

			int pointer = CIGAR_X;
			if (read_char_cache == scaff[ref_index]) {

				diag_cell += mat;

				pointer = CIGAR_EQ;
			} else if (read_char_cache != 'N' && read_char_cache != line_end) {
				diag_cell += mis;
			}

			short up_cell = local_mat_line[ref_index + 1] + gap_read;

			//find max
			short max_cell = 0;
			max_cell = max(left_cell, max_cell);
			max_cell = max(diag_cell, max_cell);
			max_cell = max(up_cell, max_cell);

			if (max_cell == up_cell
					&& max_cell != (local_mat_line[ref_index] + mis)) {
				//pointer = 2;
				pointer = CIGAR_I;
			} else if (max_cell == left_cell
					&& max_cell != (local_mat_line[ref_index] + mis)) {
				//pointer = 1;
				pointer = CIGAR_D;
			} else if (max_cell > 0
					&& (max_cell == diag_cell
							|| (max_cell == local_mat_line[ref_index] + mis
									|| max_cell
											== local_mat_line[ref_index] + mat))) {
				//pointer = 4;
			} else {
				pointer = CIGAR_STOP;
			}

			matrix[(ref_index + 1)].direction = pointer;
			matrix[(ref_index + 1)].score = max_cell;
			matrix[(ref_index + 1)].indelRun = 0;

			if (max_cell > curr_max) {
				curr_max = max_cell;
				fwResults[param_best_ref_index] = ref_index;
				fwResults[param_best_read_index] = read_index;
				fwResults[3] = curr_max;

			}
			left_cell = max_cell;
			local_mat_line[ref_index] = max_cell;
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
}

int SWCPUCor::printCigarElement(char const op, short const length,
		char * cigar) {
	int offset = 0;
	offset = sprintf(cigar, "%d%c", length, op);
	return offset;
}

void SWCPUCor::computeCigarMD(Align & result, int const gpuCigarOffset,
		short const * const gpuCigar, char const * const refSeq,
		int corr_length, int read_length) {
	int alignment_length = corr_length + read_length + 1;

	int cigar_offset = 0;
	int md_offset = 0;

	if ((gpuCigar[gpuCigarOffset] >> 4) > 0) {
		cigar_offset += printCigarElement('S', gpuCigar[gpuCigarOffset] >> 4,
				result.pRef + cigar_offset);
		result.QStart = gpuCigar[gpuCigarOffset] >> 4;
	}

	int cigar_m_length = 0;
	int md_eq_length = 0;
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
			md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
			for (int k = 0; k < length; ++k) {
				md_offset += sprintf(result.pQry + md_offset, "%c",
						refSeq[ref_index++]);
			}
			md_eq_length = 0;

			break;
		case CIGAR_EQ:
			cigar_m_length += length;
			md_eq_length += length;
			ref_index += length;
			break;
		case CIGAR_D:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length,
						result.pRef + cigar_offset);
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('D', length,
					result.pRef + cigar_offset);

			md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
			md_eq_length = 0;
			result.pQry[md_offset++] = '^';
			for (int k = 0; k < length; ++k) {
				result.pQry[md_offset++] = refSeq[ref_index++];
			}

			break;
		case CIGAR_I:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length,
						result.pRef + cigar_offset);
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('I', length,
					result.pRef + cigar_offset);

			break;
		default:
			Log.Error("Invalid cigar string: %d", op);
			std::cout << "Offset: " << gpuCigarOffset << std::endl;
			for (int x = 0; x < alignment_length * 2; ++x) {
				std::cout << gpuCigar[x] << " ";
			}
			std::cout << std::endl;
			exit(1);
		}
	}
	md_offset += sprintf(result.pQry + md_offset, "%d", md_eq_length);
	if (cigar_m_length > 0) {
		cigar_offset += printCigarElement('M', cigar_m_length,
				result.pRef + cigar_offset);
		cigar_m_length = 0;
	}

	if ((gpuCigar[alignment_length - 1] >> 4) > 0) {
		cigar_offset += printCigarElement('S',
				gpuCigar[alignment_length - 1] >> 4,
				result.pRef + cigar_offset);
		result.QEnd = gpuCigar[alignment_length - 1] >> 4;
	}

	result.pRef[cigar_offset] = '\0';
	result.pQry[md_offset] = '\0';
}

void SWCPUCor::Backtracking_CIGAR(char const * const scaff,
		char const * const read, short *& fwdResults, short *& alignments,
		int corr_length, int read_length, int alignment_length,
		MatrixElement * mat_pointer) {

	MatrixElement * matrix = mat_pointer;

	short best_read_index = fwdResults[param_best_read_index];
	short best_ref_index = fwdResults[param_best_ref_index];

	if (best_read_index > 0) {
		matrix += (((corr_length + 1) * (best_read_index + 1)));

		short abs_ref_index = best_ref_index + best_read_index;
		short alignment_index = alignment_length - 1;

		int pointer = CIGAR_STOP;
		int cigar_element = CIGAR_S;
		int cigar_length = fwdResults[qend];
		while ((pointer = matrix[(best_ref_index + 1)].direction) != CIGAR_STOP) {

			if (pointer == CIGAR_X || pointer == CIGAR_EQ) {
				matrix -= ((corr_length + 1));
				best_read_index -= 1;
				abs_ref_index -= 1;
			} else if (pointer == CIGAR_I) {
				matrix -= ((corr_length + 1));
				best_read_index -= 1;
				best_ref_index += 1;
			} else {
				best_ref_index -= 1;
				abs_ref_index -= 1;
			}

			if (pointer == cigar_element) {
				cigar_length += 1;
			} else {
				alignments[alignment_index--] = (cigar_length << 4
						| cigar_element);
				Log.Message("%d-%d", cigar_element, cigar_length);
				cigar_element = pointer;
				cigar_length = 1;
			}
		}
		alignments[alignment_index--] = (cigar_length << 4 | cigar_element);
		alignments[alignment_index] = ((best_read_index + 1) << 4 | CIGAR_S);
		fwdResults[ref_position] = abs_ref_index + 1;
		fwdResults[qstart] = best_read_index + 1;
		//qend was set by "forward" kernel
		fwdResults[alignment_offset] = alignment_index;
	}
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

	Log.Error("Unsupported alignment mode %i", mode);
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
			printf("%*d ", 2, cell->score);
		}
		printf("\n");
	}
}

int SWCPUCor::SingleAlign(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, Align & align,
		void * extData) {

	Log.Message("Aligning: ");
	Log.Message("%s", refSeq);
	Log.Message("%s", qrySeq);

	int corr_length = corridor;
	int read_length = strlen(qrySeq);
	int alignment_length = (corr_length + read_length + 1);

	int mem_matrix = (Config.GetInt("corridor") + 2)
			* (Config.GetInt("qry_max_len") + 1);

	MatrixElement * mat_pointer = new MatrixElement[mem_matrix];

	short * local_mat_line = new short[corr_length];

	memset(mat_pointer, CIGAR_STOP, mem_matrix * sizeof(char));

	short * fwdResults = new short[result_number];
	SW_Score(refSeq, qrySeq, fwdResults, corr_length, mat_pointer,
			local_mat_line);

	print_matrix(alignment_length, refSeq, read_length, qrySeq, corr_length,
			mat_pointer);
	Log.Message("%d, %d, %d, %d", fwdResults[0], fwdResults[1], fwdResults[2], fwdResults[3]);

	align.pBuffer1 = new char[read_length * 4];
	align.pBuffer2 = new char[read_length * 4];

	short * alignments = new short[alignment_length * 2];
	memset(alignments, '\0', alignment_length * 2 * sizeof(short));

//	int si = (Config.GetInt("corridor") + 2)
//			* (Config.GetInt("qry_max_len") + 1);

	Backtracking_CIGAR(refSeq, qrySeq, fwdResults, alignments, corr_length,
			read_length, alignment_length, mat_pointer);

	computeCigarMD(align, fwdResults[3], alignments, refSeq + fwdResults[0],
			corr_length, read_length);
	align.PositionOffset = fwdResults[0];

	delete[] alignments;
	delete[] fwdResults;
	delete[] mat_pointer;
	delete[] local_mat_line;

	return 1;

}

int SWCPUCor::BatchScore(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, float * const results,
		void * extData) {
	throw "Not implemented";
	return 0;

}


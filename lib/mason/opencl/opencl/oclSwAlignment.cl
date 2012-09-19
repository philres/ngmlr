#ifdef __GPU__
__kernel void oclSW_Backtracking(__global char const * scaff, __global char const * read, __global short * result, __global char * matrix, __global char * alignments) {

	matrix = matrix + ((global_index / threads_per_block) * threads_per_block * ((corridor_length + 1) * (read_length + 1))) + global_index % threads_per_block;
	read = read + (global_index * read_length);
	scaff = scaff + ((global_index / interleave_number) * interleave_number * ref_length) + global_index % interleave_number;
	result = result + result_number * global_index;

	//Index index;
	short best_read_index = result[param_best_read_index];
	short best_ref_index = result[param_best_ref_index];

	if (best_read_index > 0) {

		alignments = alignments + (global_index * alignment_length * 2);
		matrix += (((corridor_length + 1) * (best_read_index + 1)) * threads_per_block);

		short abs_ref_index = best_ref_index + best_read_index;
		short alignment_index = alignment_length - 2;

		alignments[alignment_length - 1] = '\0';
		alignments[alignment_length - 1 + alignment_length] = '\0';

		int pointer = CIGAR_STOP;
		while ((pointer = matrix[(best_ref_index + 1) * threads_per_block]) != CIGAR_STOP) {
			//printf("%d ", pointer);
			if (pointer == CIGAR_X || pointer == CIGAR_EQ) {
				alignments[alignment_index] = scaff[abs_ref_index-- * interleave_number];
				alignments[alignment_index + alignment_length] = read[best_read_index];
				matrix -= ((corridor_length + 1) * threads_per_block);
				best_read_index -= 1;
			} else if (pointer == CIGAR_I) {
				alignments[alignment_index] = '-';
				alignments[alignment_index + alignment_length] = read[best_read_index];
				matrix -= ((corridor_length + 1) * threads_per_block);
				best_read_index -= 1;
				best_ref_index += 1;
			} else {
				alignments[alignment_index] = scaff[abs_ref_index-- * interleave_number];
				alignments[alignment_index + alignment_length] = '-';
				best_ref_index -= 1;
			}
			//printf("%c %c\n", alignments[alignment_index], alignments[alignment_index + alignment_length]);
			alignment_index -= 1;
		}
		//printf("%d ", pointer);
		result[ref_position] = abs_ref_index + 1;
		result[qstart] = best_read_index + 1;
		//qend was set by "forward" kernel
		result[alignment_offset] = alignment_index;
		//printf("ref_position: %d, qstart: %d, alignment_offset: %d, alignment_length: %d", result[ref_position], result[qstart], result[alignment_offset], alignment_length);
	}

}
#endif

#ifdef __CPU__
__kernel void oclSW_Backtracking(__global char const * scaff, __global char const * read, __global short * result, __global char * _matrix, __global char * alignments) {

	//	matrix = matrix + ((global_index / threads_per_block) * threads_per_block * ((corridor_length + 1) * (read_length + 1))) + global_index % threads_per_block;

	read = read + (global_index * read_length * 4);
	scaff = scaff + (global_index * ref_length * 4);
	result = result + result_number * global_index * 4;
	alignments = alignments + (global_index * alignment_length * 2) * 4;
	_matrix = _matrix + (global_index * ((corridor_length + 1) * (read_length + 1)) * 4);

	int step_size = ((corridor_length + 1) * 4);

	for (int k = 0; k < 4; ++k) {

		__global char * matrix = _matrix;

		//Index index;
		short best_read_index = result[param_best_read_index * 4 + k];
		short best_ref_index = result[param_best_ref_index * 4 + k];

		if (best_read_index > 0) {

			//		matrix += (((corridor_length + 1) * (best_read_index + 1)) * threads_per_block);
			matrix += step_size * (best_read_index + 1);

			short abs_ref_index = best_ref_index + best_read_index;
			short alignment_index = alignment_length - 2;

			alignments[alignment_length - 1] = '\0';
			alignments[alignment_length - 1 + alignment_length] = '\0';

			int pointer = CIGAR_STOP;
			//		while ((pointer = matrix[(best_ref_index + 1) * threads_per_block]) != CIGAR_STOP) {
			while ((pointer = matrix[(best_ref_index + 1) * 4 + k]) != CIGAR_STOP) {
				//printf("%d ", pointer);
				if (pointer == CIGAR_X || pointer == CIGAR_EQ) {
					alignments[alignment_index] = scaff[abs_ref_index--];
					alignments[alignment_index + alignment_length] = read[best_read_index];
					//				matrix -= ((corridor_length + 1) * threads_per_block);
					matrix -= step_size;
					best_read_index -= 1;
				} else if (pointer == CIGAR_I) {
					alignments[alignment_index] = '-';
					alignments[alignment_index + alignment_length] = read[best_read_index];
					//				matrix -= ((corridor_length + 1) * threads_per_block);
					matrix -= step_size;
					best_read_index -= 1;
					best_ref_index += 1;
				} else {
					alignments[alignment_index] = scaff[abs_ref_index--];
					alignments[alignment_index + alignment_length] = '-';
					best_ref_index -= 1;
				}
				alignment_index -= 1;
			}
			result[ref_position * 4 + k] = abs_ref_index + 1;
			result[qstart * 4 + k] = best_read_index + 1;
			//qend was set by "forward" kernel
			result[alignment_offset * 4 + k] = alignment_index;
		}

		read = read + read_length;
		scaff = scaff + ref_length;
		//		result = result + result_number;
		alignments = alignments + alignment_length * 2;
		//		_matrix = _matrix + ((corridor_length + 1) * (read_length + 1));
	}

//	short tmp = result[1];
//	result[1] = result[4];
//	result[4] = tmp;
//
//	tmp = result[2];
//	result[2] = result[8];
//	result[8] = tmp;
//
//	tmp = result[3];
//	result[3] = result[12];
//	result[12] = tmp;
//
//	tmp = result[6];
//	result[9] = result[4];
//	result[9] = tmp;
//
//	tmp = result[7];
//	result[7] = result[13];
//	result[13] = tmp;
//
//	tmp = result[11];
//	result[11] = result[14];
//	result[14] = tmp;

}
#endif

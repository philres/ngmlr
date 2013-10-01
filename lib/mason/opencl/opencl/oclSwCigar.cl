#ifdef __GPU__
__kernel void oclSW_Backtracking(__global char const * scaff, __global char const * read, __global short * result, __global char * matrix, __global short * alignments) {

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
		short alignment_index = alignment_length - 1;

		int pointer = CIGAR_STOP;
		int cigar_element = CIGAR_S;
		int cigar_length = result[qend];
		while ((pointer = matrix[(best_ref_index + 1) * threads_per_block]) != CIGAR_STOP) {
			//printf("%d ", pointer);
			if (pointer == CIGAR_X || pointer == CIGAR_EQ) {
				matrix -= ((corridor_length + 1) * threads_per_block);
				best_read_index -= 1;
				abs_ref_index -= 1;
			} else if (pointer == CIGAR_I) {
				matrix -= ((corridor_length + 1) * threads_per_block);
				best_read_index -= 1;
				best_ref_index += 1;
			} else {
				best_ref_index -= 1;
				abs_ref_index -= 1;
			}

			if (pointer == cigar_element) {
				cigar_length += 1;
			} else {
				alignments[alignment_index--] = (cigar_length << 4 | cigar_element);
				cigar_element = pointer;
				cigar_length = 1;
			}
		}
		alignments[alignment_index--] = (cigar_length << 4 | cigar_element);
		alignments[alignment_index] = ((best_read_index + 1) << 4 | CIGAR_S);

		result[ref_position] = abs_ref_index + 1;
		result[qstart] = best_read_index + 1;
		//qend was set by "forward" kernel
		result[alignment_offset] = alignment_index;
	}

}
#endif

#ifdef __CPU__
__kernel void oclSW_Backtracking(__global char const * scaff, __global char const * read, __global short * result, __global char * _matrix, __global short * alignments) {

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

				matrix += step_size * (best_read_index + 1);

				short abs_ref_index = best_ref_index + best_read_index;
				short alignment_index = alignment_length - 1;

				int pointer = CIGAR_STOP;
				int cigar_element = CIGAR_S;
				int cigar_length = result[qend * 4 + k];
				while ((pointer = matrix[(best_ref_index + 1) * 4 + k]) != CIGAR_STOP) {
					//printf("%d ", pointer);
					if (pointer == CIGAR_X || pointer == CIGAR_EQ) {
						matrix -= step_size;
						best_read_index -= 1;
						abs_ref_index -= 1;
					} else if (pointer == CIGAR_I) {
						matrix -= step_size;
						best_read_index -= 1;
						best_ref_index += 1;
					} else {
						best_ref_index -= 1;
						abs_ref_index -= 1;
					}

					if (pointer == cigar_element) {
						cigar_length += 1;
					} else {
						alignments[alignment_index--] = (cigar_length << 4 | cigar_element);
						cigar_element = pointer;
						cigar_length = 1;
					}
				}
				alignments[alignment_index--] = (cigar_length << 4 | cigar_element);
				alignments[alignment_index] = ((best_read_index + 1) << 4 | CIGAR_S);

				result[ref_position * 4 + k] = abs_ref_index + 1;
				result[qstart * 4 + k] = best_read_index + 1;
				//qend was set by "forward" kernel
				result[alignment_offset * 4 + k] = alignment_index;
			}

		read = read + read_length;
		scaff = scaff + ref_length;
		alignments = alignments + alignment_length * 2;
	}
}
#endif


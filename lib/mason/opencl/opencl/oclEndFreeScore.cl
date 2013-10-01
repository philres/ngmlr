#ifdef __CPU__

__kernel void oclSW_Global(__global char const * scaff, __global char const * read, __global float * results) {

	scaff = scaff + (global_index * ref_length * 4);
	read = read + (global_index * read_length * 4);

	float4 l_matrix_lines[MATRIX_LENGTH];
	float4 * matrix_lines = l_matrix_lines + local_index * corridor_length;

	float4 curr_max = (float4)(short_min, short_min, short_min, short_min);

	if (read[0] != line_end) {

		//Init matrix lines
		for (short i = 0; i < corridor_length; ++i) {
			matrix_lines[i] = null4;
		}
		matrix_lines[(corridor_length - 1)] = short_min;

		for (int read_pos = 0; read_pos < read_length; ++read_pos) {

			int4 read_char_cache = (int4)(trans[read[read_pos]], trans[read[read_pos + 1 * read_length]], trans[read[read_pos + 2 * read_length]], trans[read[read_pos + 3 * read_length]]);

			float4 left_cell = short_min;
			for (int ref_index = 0; ref_index < corridor_length - 1; ++ref_index) {

				//Max can be set to left_cell because we are computing the score row-wise, so the actual max is the left_cell for the next cell
				left_cell = max(matrix_lines[(ref_index + 1)] + gap_read, left_cell + gap_ref);

				float4 scores4 = (float4)(scores[read_char_cache.s0][trans[scaff[ref_index]]], scores[read_char_cache.s1][trans[scaff[ref_index + 1 * ref_length]]], scores[read_char_cache.s2][trans[scaff[ref_index + 2 * ref_length]]], scores[read_char_cache.s3][trans[scaff[ref_index + 3 * ref_length]]]);
				left_cell = max(matrix_lines[ref_index] + scores4, left_cell);

				matrix_lines[ref_index] = left_cell;
			}
			scaff += 1;
		}

		for (short i = 0; i < corridor_length; ++i) {
			curr_max = max(curr_max, matrix_lines[i]);
		}
	}
	vstore4(curr_max, global_index, results);
}

__kernel void oclSW_ScoreGlobal(__global char const * scaff, __global char const * read, __global short * result, __global char * matrix) {

	matrix = matrix + (global_index * ((corridor_length + 1) * (read_length + 1)) * 4);
	read = read + (global_index * read_length * 4);
	scaff = scaff + (global_index * ref_length * 4);
	result = result + result_number * global_index * 4;


	float4 best_ref_index4 = null4;
	float4 read_index = null4;
	
	if (read[0] != line_end) {
		float4 l_matrix_lines[MATRIX_LENGTH];
		float4 * local_matrix_line = l_matrix_lines + local_index * corridor_length;
		//Init matrix lines
		for (short i = 0; i < corridor_length; ++i) {
			local_matrix_line[i] = null4;
			vstore4(CIGAR_STOP, i, matrix);
		}
		vstore4(CIGAR_STOP, corridor_length, matrix);
		local_matrix_line[corridor_length - 1] = short_min;		

		for (short read_pos = 0; read_pos < read_length; ++read_pos) {

			int4 read_char_cache = (int4)(trans[read[read_pos + 0 * read_length]], trans[read[read_pos + 1 * read_length]], trans[read[read_pos + 2 * read_length]], trans[read[read_pos + 3 * read_length]]);

			matrix += (corridor_length + 1) * 4;

			float4 left_cell = short_min;

			//When backtracking hits corridor border => left_up
			vstore4(CIGAR_X, 0, matrix);

			for (short ref_index = 0; ref_index < corridor_length - 1; ++ref_index) {

				//init values
				left_cell += gap_ref4;

				float4 score = (float4)(scores[read_char_cache.s0][trans[scaff[ref_index]]], scores[read_char_cache.s1][trans[scaff[ref_index + 1 * ref_length]]], scores[read_char_cache.s2][trans[scaff[ref_index + 2 * ref_length]]], scores[read_char_cache.s3][trans[scaff[ref_index + 3 * ref_length]]]);
				float4 diag_cell = local_matrix_line[ref_index] + score;
				float4 pointerX = select(CIGAR_X4, CIGAR_EQ4, (score == match));

				float4 up_cell = local_matrix_line[(ref_index + 1)] + gap_read4;

				//find max
				float4 max_cell = max(left_cell, diag_cell);
				max_cell = max(up_cell, max_cell);

				float4 pointer = 0;
				pointer = select(pointer, CIGAR_D, max_cell == left_cell);
				pointer = select(pointer, CIGAR_I, max_cell == up_cell);
				pointer = select(pointer, pointerX, max_cell == diag_cell || max_cell == (local_matrix_line[ref_index] + mismatch));

				vstore4(convert_char4(pointer), (ref_index + 1), matrix);

				left_cell = max_cell;
				local_matrix_line[ref_index] = max_cell;
			}

			vstore4(CIGAR_X, corridor_length, matrix);
			scaff += 1;
			//Find actual read length. Only increase when not \0 == 6 (see oclDefines.cl)
			read_index = select(read_index + 1, read_index, read_char_cache == (int4) 6);
		}

		float4 curr_max = short_min;
		for (short i = 0; i < corridor_length; ++i) {
			int4 lg = local_matrix_line[i] > curr_max;
			curr_max = select(curr_max, local_matrix_line[i], lg);
			best_ref_index4 = select(best_ref_index4, i, lg);
		}
	}

	vstore4(convert_short4(read_index - 1), param_best_read_index, result);
	vstore4(convert_short4(best_ref_index4), param_best_ref_index, result);
	vstore4(convert_short4(null4), qend, result);
}

#endif

#ifdef __GPU__

__kernel void oclSW_Global(__global char const * scaff, __global char const * read, __global float * results) {

	read = read + (global_index * read_length);
	scaff = scaff + ((global_index / interleave_number) * interleave_number * ref_length) + (global_index % interleave_number);

	__local
	short l_matrix_lines[MATRIX_LENGTH];
	__local
	short * matrix_lines = l_matrix_lines + local_index;
	//Init matrix lines
	for (short i = 0; i < corridor_length; ++i) {
		matrix_lines[i * threads_per_block] = 0;
	}
	matrix_lines[(corridor_length - 1) * threads_per_block] = short_min;

	for (int read_pos = 0; read_pos < read_length; ++read_pos) {
		char read_char_cache = *read++;

		short left_cell = short_min;
		for (int ref_index = 0; ref_index < (corridor_length - 1); ++ref_index) {

			short diag_cell = matrix_lines[ref_index * threads_per_block];
			diag_cell += scores[trans[read_char_cache]][trans[scaff[ref_index * interleave_number]]];

			//Max can be set to left_cell because we are computing the score row-wise, so the actual max is the left_cell for the next cell
			left_cell = max((matrix_lines[ref_index * threads_per_block + threads_per_block] + gap_read), left_cell + gap_ref);
			left_cell = max(diag_cell, left_cell);

			matrix_lines[ref_index * threads_per_block] = left_cell;
		}
		scaff += interleave_number;
	}
	short curr_max = short_min;
	for (short i = 0; i < corridor_length; ++i) {
		curr_max = max(curr_max, matrix_lines[i * threads_per_block]);
	}
	results[global_index] = curr_max;
}

__kernel void oclSW_ScoreGlobal(__global char const * scaff, __global char const * read, __global short * result, __global char * matrix) {

	matrix = matrix + ((global_index / threads_per_block) * threads_per_block * ((corridor_length + 1) * (read_length + 1))) + global_index % threads_per_block;
	read = read + (global_index * read_length);
	scaff = scaff + ((global_index / interleave_number) * interleave_number * ref_length) + global_index % interleave_number;
	result = result + result_number * global_index;

	__local
	short l_matrix_lines[MATRIX_LENGTH];
	__local
	short * local_matrix_line = l_matrix_lines + local_index;

	//Init matrix lines
	for (short i = 0; i < corridor_length; ++i) {
		local_matrix_line[i * threads_per_block] = 0;
		matrix[i * threads_per_block] = CIGAR_STOP;
	}
	matrix[corridor_length * threads_per_block] = CIGAR_STOP;

	local_matrix_line[(corridor_length - 1)] = short_min;

	short read_index = 0;

//			bool test = read[0] != '\0';	
//			if (test) {
//				//printf("GPU Ref : %s\n", scaff);
//				//printf("GPU Read: %s\n", read);
//				printf("S:\t");
//				char c;
//				int i = 0;
//				while ((c = scaff[i * interleave_number]) != '\0') {
//					printf("\t%c", c);
//					i += 1;
//				}
//				printf("\n");
//			}


	for (char read_char_cache; (read_char_cache = *read) != line_end; ++read) {
//										printf("%c:\t", read_char_cache);
//										for (short i = 0; i < read_index; ++i) {
//											printf("\t");
//										}
		matrix += (corridor_length + 1) * threads_per_block;

		short left_cell = short_min;
		matrix[0] = CIGAR_X;
		for (short ref_index = 0; ref_index < corridor_length - 1; ++ref_index) {

			//init values
			left_cell += gap_ref;
			short diag_cell = local_matrix_line[ref_index * threads_per_block];
			//			printf("%c == %c\n", read_char_cache, scaff[ref_index * interleave_number]);
			int pointer = CIGAR_X;

			diag_cell += scores[trans[read_char_cache]][trans[scaff[ref_index * interleave_number]]];
			pointer = select(CIGAR_X, CIGAR_EQ, (read_char_cache == scaff[ref_index * interleave_number]));

			short up_cell = local_matrix_line[(ref_index + 1) * threads_per_block] + gap_read;

			//find max
			short max_cell = 0;

			max_cell = max(diag_cell, left_cell);
			max_cell = max(up_cell, max_cell);


			if(max_cell == diag_cell ||
				 	   max_cell == (local_matrix_line[ref_index * threads_per_block] + mismatch)) {
			} else if(max_cell == up_cell) {
				pointer = CIGAR_I;
			} else {
				pointer = CIGAR_D;
			}

//			printf("\t%d", (int)max_cell);
			matrix[(ref_index + 1) * threads_per_block] = pointer;

			left_cell = max_cell;
			local_matrix_line[ref_index * threads_per_block] = max_cell;
		}
		matrix[corridor_length * threads_per_block] = CIGAR_X;
		scaff += interleave_number;
		read_index += 1;

//				for (short i = 0; i < corridor_length + 1; ++i) {
//					printf("\t%d ", matrix[i* threads_per_block]);
//				}
		//printf("\n");
	}

	short curr_max = short_min;
	for (short i = 0; i < corridor_length; ++i) {
		//		printf("\t%d ", local_matrix_line[i]);
		if (local_matrix_line[i * threads_per_block] > curr_max) {
			curr_max = local_matrix_line[i * threads_per_block];
			result[param_best_read_index] = read_index - 1;
			result[param_best_ref_index] = i;
		}
	}

	result[qend] = 0;
	if (read_index == 0) {
		result[param_best_read_index] = result[param_best_read_index] = 0;
	}
}

#endif

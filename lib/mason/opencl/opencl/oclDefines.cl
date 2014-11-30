#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store: enable

#ifdef __NVIDIA__
#define __cl_constant __constant const
#else
#define __cl_constant __constant
#endif

#define global_index get_global_id(0)
#define local_index get_local_id(0)

#define ref_position 0
#define qstart 1
#define qend 2
#define alignment_offset 3

#define param_best_read_index 0
#define param_best_ref_index 1

#define CIGAR_STOP 10


__cl_constant char line_end = '\0';
__cl_constant short null = 0;

//Use FLT_MIN instead
#define short_min -16000

#ifdef __CPU__

__cl_constant float4 match4 = (float4)(match);
__cl_constant float4 mismatch4 = (float4)(mismatch);

__cl_constant float4 gap_read4 = (float4)(gap_read);
__cl_constant float4 gap_ref4 = (float4)(gap_ref);

__cl_constant float4 null4 = (float4)(0);

__cl_constant float4 CIGAR_X4 = (float4)(CIGAR_X);
__cl_constant float4 CIGAR_EQ4 = (float4)(CIGAR_EQ);
__cl_constant float4 CIGAR_I4 = (float4)(CIGAR_I);
__cl_constant float4 CIGAR_D4 = (float4)(CIGAR_D);
__cl_constant float4 CIGAR_STOP4 = (float4)(CIGAR_STOP);



__cl_constant short4 match4s = (short4)(match);
__cl_constant short4 mismatch4s = (short4)(mismatch);

__cl_constant short4 gap_read4s = (short4)(gap_read);
__cl_constant short4 gap_ref4s = (short4)(gap_ref);

__cl_constant short4 null4s = (short4)(0);

__cl_constant short4 CIGAR_X4s = (short4)(CIGAR_X);
__cl_constant short4 CIGAR_EQ4s = (short4)(CIGAR_EQ);
__cl_constant short4 CIGAR_I4s = (short4)(CIGAR_I);
__cl_constant short4 CIGAR_D4s = (short4)(CIGAR_D);
__cl_constant short4 CIGAR_STOP4s = (short4)(CIGAR_STOP);

#endif

__cl_constant int trans[256] =
		{ 		6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 5, 4,
				4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 5, 4,
				4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
				4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };


									//Ref ->
									//A			C		  G         T         X         N         \0
__cl_constant float scores[7][7] = { {match,    mismatch, mismatch, mismatch, mismatch, mismatch, 0},
/*Read*/							 {mismatch, match,    mismatch, mismatch, mismatch, mismatch, 0},
/*|*/								 {mismatch, mismatch, match,    mismatch, mismatch, mismatch, 0},
/*v*/								 {mismatch, mismatch, mismatch, match,    mismatch, mismatch, 0},
									 {mismatch, mismatch, mismatch, mismatch, mismatch, mismatch, 0},
									 {0,        0,        0,        0,        mismatch, mismatch, mismatch},
									 {0,        0,        0,        0,        0,        0       , 0} };

									   //Ref ->
									   //A			C		  G         T         X         N         \0
__cl_constant float scoresTC[7][7] = { {match,    mismatch,   mismatch, mismatch, mismatch, mismatch, 0},
/*Read*/							   {mismatch, match,      mismatch, mismatch, mismatch, mismatch, 0},
/*|*/								   {mismatch, mismatch,   match,    mismatch, mismatch, mismatch, 0},
/*v*/								   {mismatch, mismatchBS, mismatch, matchBS,  mismatch, mismatch, 0},
 	 	 	 	 	 	 	 	 	   {mismatch, mismatch,   mismatch, mismatch, mismatch, mismatch, 0},
 	 	 	 	 	 	 	 	 	   {0,        0,          0,        0,        mismatch, mismatch, mismatch},
 	 	 	 	 	 	 	 	 	   {0,        0,          0,        0,        0,        0       , 0} };


__cl_constant float scoresAG[7][7] = { {matchBS,  mismatch, mismatchBS, mismatch, mismatch, mismatch, 0},
/*Read*/							   {mismatch, match,    mismatch,   mismatch, mismatch, mismatch, 0},
/*|*/								   {mismatch, mismatch, match,      mismatch, mismatch, mismatch, 0},
/*v*/								   {mismatch, mismatch, mismatch,   match,    mismatch, mismatch, 0},
 	 	 	 	 	 	 	 	 	   {mismatch, mismatch, mismatch,   mismatch, mismatch, mismatch, 0},
 	 	 	 	 	 	 	 	 	   {0,        0,        0,          0,        0,        0       , 0}       ,
 	 	 	 	 	 	 	 	 	   {0,        0,        0,          0,        0,        0       , 0} };

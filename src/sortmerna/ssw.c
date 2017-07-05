/* The MIT License

   Copyright (c) 2012-1015 Boston College.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.	
*/

/* Contact: Mengyao Zhao <zhangmp@bc.edu> */

/*
 *  ssw.c
 *
 *  Created by Mengyao Zhao on 6/22/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *	Version 0.1.4
 *	Last revision by Mengyao Zhao on 12/07/12.
 *
 */

/*
 * Important: parts of this file have been modified by Evguenia Kopylova
 *            to work with SortMeRNA. Please download the code published
 *            with the SSW library for the original version.
 *            July 10, 2013
 */

#include <emmintrin.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../include/ssw.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/* Convert the coordinate in the scoring matrix into the coordinate in one line of the band. */
#define set_u(u, w, i, j) { int x=(i)-(w); x=x>0?x:0; (u)=(j)-x+1; }

/* Convert the coordinate in the direction matrix into the coordinate in one line of the band. */
#define set_d(u, w, i, j, p) { int x=(i)-(w); x=x>0?x:0; x=(j)-x; (u)=x*3+p; }

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

typedef struct {
	uint16_t score;
	int32_t ref;	 //0-based position 
	int32_t read;    //alignment ending position on read, 0-based 
} alignment_end;

typedef struct {
	uint32_t* seq;
	int32_t length;
} cigar;

struct _profile{
	__m128i* profile_byte;	// 0: none
	__m128i* profile_word;	// 0: none
	const int8_t* read;
	const int8_t* mat;
	int32_t readLen;
	int32_t n;
	uint8_t bias;
};

int8_t rc_table[128] = {
    4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 4, 4,
    4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 4, 4,
    4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch. */
__m128i* qP_byte (const int8_t* read_num,
				  const int8_t* mat,
				  const int32_t readLen,
				  const int32_t n,	/* the edge length of the square matrix mat */
				  uint8_t bias) {
 
	int32_t segLen = (readLen + 15) / 16; /* Split the 128 bit register into 16 pieces. 
								     Each piece is 8 bit. Split the read into 16 segments. 
								     Calculate 16 segments in parallel.
								   */
	__m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
	int8_t* t = (int8_t*)vProfile;
	int32_t nt, i, j, segNum;
	
	/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch */
	for (nt = 0; LIKELY(nt < n); nt ++) {
		for (i = 0; i < segLen; i ++) {
			j = i; 
			for (segNum = 0; LIKELY(segNum < 16) ; segNum ++) {
				*t++ = j>= readLen ? bias : mat[nt * n + read_num[j]] + bias;
				j += segLen;
			}
		}
	}
	return vProfile;
}

/* Striped Smith-Waterman
   Record the highest score of each reference position. 
   Return the alignment score and ending position of the best alignment, 2nd best alignment, etc. 
   Gap begin and gap extension are different. 
   wight_match > 0, all other weights < 0.
   The returned positions are 0-based.
 */ 
alignment_end* sw_sse2_byte (const int8_t* ref,
							 int8_t ref_dir,	// 0: forward ref; 1: reverse ref
							 int32_t refLen,
							 int32_t readLen, 
							 const uint8_t weight_gapO, /* will be used as - */
							 const uint8_t weight_gapE, /* will be used as - */
							 __m128i* vProfile,
							 uint8_t terminate,	/* the best alignment score: used to terminate 
												   the matrix calculation when locating the 
												   alignment beginning point. If this score 
												   is set to 0, it will not be used */
	 						 uint8_t bias,  /* Shift 0 point to a positive value. */
							 int32_t maskLen) {  
      
#define max16(m, vm) (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 8)); \
					  (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 4)); \
					  (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 2)); \
					  (vm) = _mm_max_epu8((vm), _mm_srli_si128((vm), 1)); \
					  (m) = _mm_extract_epi16((vm), 0)

	uint8_t max = 0;		                     /* the max alignment score */
	int32_t end_read = readLen - 1;
	int32_t end_ref = -1; /* 0_based best alignment ending point; Initialized as isn't aligned -1. */
	int32_t segLen = (readLen + 15) / 16; /* number of segment */
	
	/* array to record the largest score of each reference position */
	uint8_t* maxColumn = (uint8_t*) calloc(refLen, 1); 
	
	/* array to record the alignment read ending position of the largest score of each reference position */
	int32_t* end_read_column = (int32_t*) calloc(refLen, sizeof(int32_t));
	
	/* Define 16 byte 0 vector. */
	__m128i vZero = _mm_set1_epi32(0);

	__m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));

	int32_t i, j;
	/* 16 byte insertion begin vector */
	__m128i vGapO = _mm_set1_epi8(weight_gapO);
	
	/* 16 byte insertion extension vector */
	__m128i vGapE = _mm_set1_epi8(weight_gapE);	
	
	/* 16 byte bias vector */
	__m128i vBias = _mm_set1_epi8(bias);	

	__m128i vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
	__m128i vMaxMark = vZero; /* Trace the highest score till the previous column. */	
	__m128i vTemp;
	int32_t edge, begin = 0, end = refLen, step = 1; 
//	int32_t distance = readLen * 2 / 3;
//	int32_t distance = readLen / 2;
//	int32_t distance = readLen;

	/* outer loop to process the reference sequence */
	if (ref_dir == 1) {
		begin = refLen - 1;
		end = -1;
		step = -1;
	}
	for (i = begin; LIKELY(i != end); i += step) {
		int32_t cmp;
		__m128i e = vZero, vF = vZero, vMaxColumn = vZero; /* Initialize F value to 0. 
							   Any errors to vH values will be corrected in the Lazy_F loop. 
							 */
//		max16(maxColumn[i], vMaxColumn);
//		fprintf(stderr, "middle[%d]: %d\n", i, maxColumn[i]);

		__m128i vH = pvHStore[segLen - 1];
		vH = _mm_slli_si128 (vH, 1); /* Shift the 128-bit value in vH left by 1 byte. */
		__m128i* vP = vProfile + ref[i] * segLen; /* Right part of the vProfile */

		/* Swap the 2 H buffers. */
		__m128i* pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;
		
		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); ++j) {
			vH = _mm_adds_epu8(vH, _mm_load_si128(vP + j));
			vH = _mm_subs_epu8(vH, vBias); /* vH will be always > 0 */
	//	max16(maxColumn[i], vH);
	//	fprintf(stderr, "H[%d]: %d\n", i, maxColumn[i]);
//	int8_t* t;
//	int32_t ti;
//for (t = (int8_t*)&vH, ti = 0; ti < 16; ++ti) fprintf(stderr, "%d\t", *t++);

			/* Get max from vH, vE and vF. */
			e = _mm_load_si128(pvE + j);
			vH = _mm_max_epu8(vH, e);
			vH = _mm_max_epu8(vH, vF);
			vMaxColumn = _mm_max_epu8(vMaxColumn, vH);
			
	//	max16(maxColumn[i], vMaxColumn);
	//	fprintf(stderr, "middle[%d]: %d\n", i, maxColumn[i]);
//	for (t = (int8_t*)&vMaxColumn, ti = 0; ti < 16; ++ti) fprintf(stderr, "%d\t", *t++);

			/* Save vH values. */
			_mm_store_si128(pvHStore + j, vH);

			/* Update vE value. */
			vH = _mm_subs_epu8(vH, vGapO); /* saturation arithmetic, result >= 0 */
			e = _mm_subs_epu8(e, vGapE);
			e = _mm_max_epu8(e, vH);
			_mm_store_si128(pvE + j, e);
			
			/* Update vF value. */
			vF = _mm_subs_epu8(vF, vGapE);
			vF = _mm_max_epu8(vF, vH);
			
			/* Load the next vH. */
			vH = _mm_load_si128(pvHLoad + j);
		}

		/* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
        /* reset pointers to the start of the saved data */
        j = 0;
        vH = _mm_load_si128 (pvHStore + j);

        /*  the computed vF value is for the given column.  since */
        /*  we are at the end, we need to shift the vF value over */
        /*  to the next column. */
        vF = _mm_slli_si128 (vF, 1);
        vTemp = _mm_subs_epu8 (vH, vGapO);
		vTemp = _mm_subs_epu8 (vF, vTemp);
		vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
		cmp  = _mm_movemask_epi8 (vTemp);

        while (cmp != 0xffff) 
        {
            vH = _mm_max_epu8 (vH, vF);
			vMaxColumn = _mm_max_epu8(vMaxColumn, vH);
            _mm_store_si128 (pvHStore + j, vH);
            vF = _mm_subs_epu8 (vF, vGapE);
            j++;
            if (j >= segLen)
            {
                j = 0;
                vF = _mm_slli_si128 (vF, 1);
            }
            vH = _mm_load_si128 (pvHStore + j);

            vTemp = _mm_subs_epu8 (vH, vGapO);
            vTemp = _mm_subs_epu8 (vF, vTemp);
            vTemp = _mm_cmpeq_epi8 (vTemp, vZero);
            cmp  = _mm_movemask_epi8 (vTemp);
        }

		vMaxScore = _mm_max_epu8(vMaxScore, vMaxColumn);
		vTemp = _mm_cmpeq_epi8(vMaxMark, vMaxScore);
		cmp = _mm_movemask_epi8(vTemp);
		if (cmp != 0xffff) {
			uint8_t temp; 
			vMaxMark = vMaxScore;
			max16(temp, vMaxScore);
			vMaxScore = vMaxMark;
			
			if (LIKELY(temp > max)) {
				max = temp;
				if (max + bias >= 255) break; //overflow

				end_ref = i;
			
				/* Store the column with the highest alignment score in order to trace the alignment ending position on read. */
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}

		/* Record the max score of current column. */	
		max16(maxColumn[i], vMaxColumn);
//		fprintf(stderr, "maxColumn[%d]: %d\n", i, maxColumn[i]);
		if (maxColumn[i] == terminate) break;
	}
	
	/* Trace the alignment ending position on read. */
	uint8_t *t = (uint8_t*)pvHmax;
	int32_t column_len = segLen * 16;
	for (i = 0; LIKELY(i < column_len); ++i, ++t) {
		int32_t temp;
		if (*t == max) {
			temp = i / 16 + i % 16 * segLen;
			if (temp < end_read) end_read = temp;
		}
	}

	free(pvHmax);
	free(pvE);
	free(pvHLoad);
	free(pvHStore); 	

	/* Find the most possible 2nd best alignment. */
	alignment_end* bests = (alignment_end*) calloc(2, sizeof(alignment_end));
	bests[0].score = max + bias >= 255 ? 255 : max;
	bests[0].ref = end_ref;
	bests[0].read = end_read;
	
	bests[1].score = 0;
	bests[1].ref = 0;
	bests[1].read = 0;

	edge = (end_ref - maskLen) > 0 ? (end_ref - maskLen) : 0;
	for (i = 0; i < edge; i ++) {
//			fprintf (stderr, "maxColumn[%d]: %d\n", i, maxColumn[i]); 
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i;
		}
	}
	edge = (end_ref + maskLen) > refLen ? refLen : (end_ref + maskLen);
	for (i = edge + 1; i < refLen; i ++) {
//			fprintf (stderr, "refLen: %d\tmaxColumn[%d]: %d\n", refLen, i, maxColumn[i]);
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i;
		}
	}
    
	free(maxColumn);
	free(end_read_column);
	return bests;
}

__m128i* qP_word (const int8_t* read_num,
				  const int8_t* mat,
				  const int32_t readLen,
				  const int32_t n) { 
					
	int32_t segLen = (readLen + 7) / 8; 
	__m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
	int16_t* t = (int16_t*)vProfile;
	int32_t nt, i, j;
	int32_t segNum;
	
	/* Generate query profile rearrange query sequence & calculate the weight of match/mismatch */
	for (nt = 0; LIKELY(nt < n); nt ++) {
		for (i = 0; i < segLen; i ++) {
			j = i; 
			for (segNum = 0; LIKELY(segNum < 8) ; segNum ++) {
				*t++ = j>= readLen ? 0 : mat[nt * n + read_num[j]];
				j += segLen;
			}
		}
	}
	return vProfile;
}

alignment_end* sw_sse2_word (const int8_t* ref, 
							 int8_t ref_dir,	// 0: forward ref; 1: reverse ref
							 int32_t refLen,
							 int32_t readLen, 
							 const uint8_t weight_gapO, /* will be used as - */
							 const uint8_t weight_gapE, /* will be used as - */
						     __m128i* vProfile,
							 uint16_t terminate, 
							 int32_t maskLen) { 

#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
					(vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
					(m) = _mm_extract_epi16((vm), 0)
	
	uint16_t max = 0;		                     /* the max alignment score */
	int32_t end_read = readLen - 1;
	int32_t end_ref = 0; /* 1_based best alignment ending point; Initialized as isn't aligned - 0. */
	int32_t segLen = (readLen + 7) / 8; /* number of segment */
	
	/* array to record the largest score of each reference position */
	uint16_t* maxColumn = (uint16_t*) calloc(refLen, 2); 
	
	/* array to record the alignment read ending position of the largest score of each reference position */
	int32_t* end_read_column = (int32_t*) calloc(refLen, sizeof(int32_t));
	
	/* Define 16 byte 0 vector. */
	__m128i vZero = _mm_set1_epi32(0);

	__m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));
	__m128i* pvHmax = (__m128i*) calloc(segLen, sizeof(__m128i));

	int32_t i, j, k;
	/* 16 byte insertion begin vector */
	__m128i vGapO = _mm_set1_epi16(weight_gapO);
	
	/* 16 byte insertion extension vector */
	__m128i vGapE = _mm_set1_epi16(weight_gapE);	

	/* 16 byte bias vector */
	__m128i vMaxScore = vZero; /* Trace the highest score of the whole SW matrix. */
	__m128i vMaxMark = vZero; /* Trace the highest score till the previous column. */	
	__m128i vTemp;
	int32_t edge, begin = 0, end = refLen, step = 1;

	/* outer loop to process the reference sequence */
	if (ref_dir == 1) {
		begin = refLen - 1;
		end = -1;
		step = -1;
	}
	for (i = begin; LIKELY(i != end); i += step) {
		int32_t cmp;
		__m128i e = vZero, vF = vZero; /* Initialize F value to 0. 
							   Any errors to vH values will be corrected in the Lazy_F loop. 
							 */
		__m128i vH = pvHStore[segLen - 1];
		vH = _mm_slli_si128 (vH, 2); /* Shift the 128-bit value in vH left by 2 byte. */
		
		/* Swap the 2 H buffers. */
		__m128i* pv = pvHLoad;
		
		__m128i vMaxColumn = vZero; /* vMaxColumn is used to record the max values of column i. */
		
		__m128i* vP = vProfile + ref[i] * segLen; /* Right part of the vProfile */
		pvHLoad = pvHStore;
		pvHStore = pv;
		
		/* inner loop to process the query sequence */
		for (j = 0; LIKELY(j < segLen); j ++) {
			vH = _mm_adds_epi16(vH, _mm_load_si128(vP + j));

			/* Get max from vH, vE and vF. */
			e = _mm_load_si128(pvE + j);
			vH = _mm_max_epi16(vH, e);
			vH = _mm_max_epi16(vH, vF);
			vMaxColumn = _mm_max_epi16(vMaxColumn, vH);
			
			/* Save vH values. */
			_mm_store_si128(pvHStore + j, vH);

			/* Update vE value. */
			vH = _mm_subs_epu16(vH, vGapO); /* saturation arithmetic, result >= 0 */
			e = _mm_subs_epu16(e, vGapE);
			e = _mm_max_epi16(e, vH);
			_mm_store_si128(pvE + j, e);

			/* Update vF value. */
			vF = _mm_subs_epu16(vF, vGapE);
			vF = _mm_max_epi16(vF, vH);
			
			/* Load the next vH. */
			vH = _mm_load_si128(pvHLoad + j);
		}

		/* Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3 */
		for (k = 0; LIKELY(k < 8); ++k) {
			vF = _mm_slli_si128 (vF, 2);
			for (j = 0; LIKELY(j < segLen); ++j) {
				vH = _mm_load_si128(pvHStore + j);
				vH = _mm_max_epi16(vH, vF);
				_mm_store_si128(pvHStore + j, vH);
				vH = _mm_subs_epu16(vH, vGapO);
				vF = _mm_subs_epu16(vF, vGapE);
				if (UNLIKELY(! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH)))) goto end;
			}
		}

end:	
		vMaxScore = _mm_max_epi16(vMaxScore, vMaxColumn);	
		vTemp = _mm_cmpeq_epi16(vMaxMark, vMaxScore);
		cmp = _mm_movemask_epi8(vTemp);
		if (cmp != 0xffff) {
			uint16_t temp; 
			vMaxMark = vMaxScore;
			max8(temp, vMaxScore);
			vMaxScore = vMaxMark;
			
			if (LIKELY(temp > max)) {
				max = temp;
				end_ref = i;
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}
		
		/* Record the max score of current column. */	
		max8(maxColumn[i], vMaxColumn);
		if (maxColumn[i] == terminate) break;
	} 	

	/* Trace the alignment ending position on read. */
	uint16_t *t = (uint16_t*)pvHmax;
	int32_t column_len = segLen * 8;
	for (i = 0; LIKELY(i < column_len); ++i, ++t) {
		int32_t temp;
		if (*t == max) {
			temp = i / 8 + i % 8 * segLen;
			if (temp < end_read) end_read = temp;
		}
	}

	free(pvHmax);
	free(pvE);
	free(pvHLoad);
	free(pvHStore); 
	
	/* Find the most possible 2nd best alignment. */
	alignment_end* bests = (alignment_end*) calloc(2, sizeof(alignment_end));
	bests[0].score = max;
	bests[0].ref = end_ref;
	bests[0].read = end_read;
	
	bests[1].score = 0;
	bests[1].ref = 0;
	bests[1].read = 0;

	edge = (end_ref - maskLen) > 0 ? (end_ref - maskLen) : 0;
	for (i = 0; i < edge; i ++) {
		if (maxColumn[i] > bests[1].score) { 
			bests[1].score = maxColumn[i];
			bests[1].ref = i;
		}
	}
	edge = (end_ref + maskLen) > refLen ? refLen : (end_ref + maskLen);
	for (i = edge; i < refLen; i ++) {
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i;
		}
	}
	
	free(maxColumn);
	free(end_read_column);
	return bests;
}

cigar* banded_sw (const int8_t* ref,
				 const int8_t* read, 
				 int32_t refLen, 
				 int32_t readLen,
				 int32_t score,
				 const uint32_t weight_gapO,  /* will be used as - */
				 const uint32_t weight_gapE,  /* will be used as - */
				 int32_t band_width,
				 const int8_t* mat,	/* pointer to the weight matrix */
				 int32_t n) {
    
	uint32_t *c = (uint32_t*)malloc(16 * sizeof(uint32_t)), *c1;
	int32_t i, j, e, f, temp1, temp2, s = 16, s1 = 8, s2 = 1024, l, max = 0;
	int32_t width, width_d, *h_b, *e_b, *h_c;
	int8_t *direction, *direction_line;
	cigar* result = (cigar*)malloc(sizeof(cigar));
	h_b = (int32_t*)malloc(s1 * sizeof(int32_t)); 
	e_b = (int32_t*)malloc(s1 * sizeof(int32_t)); 
	h_c = (int32_t*)malloc(s1 * sizeof(int32_t)); 
	direction = (int8_t*)malloc(s2 * sizeof(int8_t));

	do {
		width = band_width * 2 + 3, width_d = band_width * 2 + 1;
		while (width >= s1) {
			++s1;
			kroundup32(s1);
			h_b = (int32_t*)realloc(h_b, s1 * sizeof(int32_t)); 
			e_b = (int32_t*)realloc(e_b, s1 * sizeof(int32_t)); 
			h_c = (int32_t*)realloc(h_c, s1 * sizeof(int32_t)); 
		}
		while (width_d * readLen * 3 >= s2) {
			++s2;
			kroundup32(s2);
			if (s2 < 0) {
				fprintf(stderr, "Alignment score and position are not consensus.\n");
				exit(1);
				free(h_b); //jenya
				h_b = NULL;
				free(e_b); //jenya
				e_b = NULL;
				free(h_c); //jenya
				h_c = NULL;
				free(direction); //jenya
        direction = NULL;
				free(c); //jenya
        c = NULL;
				free(result); //jenya
				return 0;
			}
			direction = (int8_t*)realloc(direction, s2 * sizeof(int8_t)); 
		}
		direction_line = direction;
		for (j = 1; LIKELY(j < width - 1); j ++) h_b[j] = 0;
		for (i = 0; LIKELY(i < readLen); i ++) {
			int32_t beg = 0, end = refLen - 1, u = 0, edge;
			j = i - band_width;	beg = beg > j ? beg : j; // band start
			j = i + band_width; end = end < j ? end : j; // band end
			edge = end + 1 < width - 1 ? end + 1 : width - 1;
			f = h_b[0] = e_b[0] = h_b[edge] = e_b[edge] = h_c[0] = 0;
			direction_line = direction + width_d * i * 3;

			for (j = beg; LIKELY(j <= end); j ++) {
				int32_t b, e1, f1, d, de, df, dh;
				set_u(u, band_width, i, j);	set_u(e, band_width, i - 1, j); 
				set_u(b, band_width, i, j - 1); set_u(d, band_width, i - 1, j - 1);
				set_d(de, band_width, i, j, 0);
				set_d(df, band_width, i, j, 1);
				set_d(dh, band_width, i, j, 2);

				temp1 = i == 0 ? -weight_gapO : h_b[e] - weight_gapO;
				temp2 = i == 0 ? -weight_gapE : e_b[e] - weight_gapE;
				e_b[u] = temp1 > temp2 ? temp1 : temp2;
				direction_line[de] = temp1 > temp2 ? 3 : 2;
		
				temp1 = h_c[b] - weight_gapO;
				temp2 = f - weight_gapE;
				f = temp1 > temp2 ? temp1 : temp2;
				direction_line[df] = temp1 > temp2 ? 5 : 4;
				
				e1 = e_b[u] > 0 ? e_b[u] : 0;
				f1 = f > 0 ? f : 0;
				temp1 = e1 > f1 ? e1 : f1;
				temp2 = h_b[d] + mat[ref[j] * n + read[i]];
				h_c[u] = temp1 > temp2 ? temp1 : temp2;
		
				if (h_c[u] > max) max = h_c[u];
		
				if (temp1 <= temp2) direction_line[dh] = 1;
				else direction_line[dh] = e1 > f1 ? direction_line[de] : direction_line[df];
			}
			for (j = 1; j <= u; j ++) h_b[j] = h_c[j];
		}
		band_width *= 2;
	} while (LIKELY(max < score));
	band_width /= 2;

	// trace back
	i = readLen - 1;
	j = refLen - 1;
	e = 0;	// Count the number of M, D or I.
	l = 0;	// record length of current cigar
	f = max = 0; // M
	temp2 = 2;	// h
	while (LIKELY(i > 0)) {
		set_d(temp1, band_width, i, j, temp2);
		switch (direction_line[temp1]) {
			case 1: 
				--i;
				--j;
				temp2 = 2;
				direction_line -= width_d * 3;
				f = 0;	// M
				break;
			case 2:
			 	--i;
				temp2 = 0;	// e
				direction_line -= width_d * 3;
				f = 1;	// I
				break;		
			case 3:
				--i;
				temp2 = 2;
				direction_line -= width_d * 3;
				f = 1;	// I
				break;
			case 4:
				--j;
				temp2 = 1;
				f = 2;	// D
				break;
			case 5:
				--j;
				temp2 = 2;
				f = 2;	// D
				break;
			default: 
				fprintf(stderr, "Trace back error: %d.\n", direction_line[temp1 - 1]); exit(1);
				return 0;
		}
		if (f == max) ++e;
		else {
			++l;
			while (l >= s) {
				++s;
				kroundup32(s);
				c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
			}
			c[l - 1] = e<<4|max;
			max = f;
			e = 1;
		}
	}
	if (f == 0) {
		++l;
		while (l >= s) {
			++s;
			kroundup32(s);
			c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
		}
		c[l - 1] = (e+1)<<4;
	}else {
		l += 2;
		while (l >= s) {
			++s;
			kroundup32(s);
			c = (uint32_t*)realloc(c, s * sizeof(uint32_t));
		}
		c[l - 2] = e<<4|f;
		c[l - 1] = 16;	// 1M
	}

	// reverse cigar
	c1 = (uint32_t*)malloc(l * sizeof(uint32_t));
	s = 0;
	e = l - 1;
	while (LIKELY(s <= e)) {			
		c1[s] = c[e];		
		c1[e] = c[s];		
		++ s;					
		-- e;						
	}								
	result->seq = c1;
	result->length = l;

	free(direction);
	direction = NULL;
	free(h_c);
	h_c = NULL;
	free(e_b);
    e_b = NULL;
    free(h_b);
    h_b = NULL;
	free(c); 
    c = NULL;
	return result;
}

int8_t* seq_reverse(const int8_t* seq, int32_t end)	/* end is 0-based alignment ending position */	
{									
	int8_t* reverse = (int8_t*)calloc(end + 1, sizeof(int8_t));	
	int32_t start = 0;
	while (LIKELY(start <= end)) {			
		reverse[start] = seq[end];		
		reverse[end] = seq[start];		
		++ start;					
		-- end;						
	}								
	return reverse;					
}
		
s_profile* ssw_init (const int8_t* read, const int32_t readLen, const int8_t* mat, const int32_t n, const int8_t score_size) {
	s_profile* p = (s_profile*)calloc(1, sizeof(struct _profile));
    if ( p == NULL )
    {
        fprintf(stderr, "    %sERROR%s: could not allocate memory for profile (ssw.c)\n","\033[0;31m","\033[0m");
        exit(EXIT_FAILURE);
    }
	p->profile_byte = 0;
	p->profile_word = 0;
	p->bias = 0;
	
	if (score_size == 0 || score_size == 2) {
		/* Find the bias to use in the substitution matrix */
		int32_t bias = 0, i;
		for (i = 0; i < n*n; i++) if (mat[i] < bias) bias = mat[i];
		bias = abs(bias);

		p->bias = bias;
		p->profile_byte = qP_byte (read, mat, readLen, n, bias);
	}
	if (score_size == 1 || score_size == 2) p->profile_word = qP_word (read, mat, readLen, n);
	p->read = read;
	p->mat = mat;
	p->readLen = readLen;
	p->n = n;
	return p;
}

void init_destroy (s_profile** p) {
    if ((*p)->profile_byte != NULL )
    {
        free((*p)->profile_byte);
        (*p)->profile_byte = NULL;
    }
    if ( (*p)->profile_word != NULL )
    {
        free((*p)->profile_word);
        (*p)->profile_word = NULL;
    }
    if (*p != NULL)
    {
        free(*p);
        *p = NULL;
    }
}

s_align* ssw_align (const s_profile* prof, 
					const int8_t* ref, 
				  	int32_t refLen, 
				  	const uint8_t weight_gapO, 
				  	const uint8_t weight_gapE, 
					const uint8_t flag,	//  (from high to low) bit 5: return the best alignment beginning position; 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar; 7: if max score >= filters, return cigar; 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
					const uint16_t filters,
					const int32_t filterd,
					const int32_t maskLen) {
    
	alignment_end* bests = 0, *bests_reverse = 0;
	__m128i* vP = 0;
	int32_t word = 0, band_width = 0, readLen = prof->readLen;
	int8_t* read_reverse = 0;
	cigar* path;
	s_align* r = (s_align*)calloc(1, sizeof(s_align));
	r->ref_begin1 = -1;
	r->read_begin1 = -1;
	r->cigar = 0;
	r->cigarLen = 0;
	if (maskLen < 15) {
		//fprintf(stderr, "When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information.\n");
	}


	// Find the alignment scores and ending positions
	if (prof->profile_byte) {

        bests = sw_sse2_byte(ref, 0, refLen, readLen, weight_gapO, weight_gapE, prof->profile_byte, -1, prof->bias, maskLen);
        
		if (prof->profile_word && bests[0].score == 255) {
			free(bests);
			bests = NULL;
			bests = sw_sse2_word(ref, 0, refLen, readLen, weight_gapO, weight_gapE, prof->profile_word, -1, maskLen);
			word = 1;
		} else if (bests[0].score == 255) {
			fprintf(stderr, "Please set 2 to the score_size parameter of the function ssw_init, otherwise the alignment results will be incorrect.\n");
			return 0;
		}
	}else if (prof->profile_word) {
		bests = sw_sse2_word(ref, 0, refLen, readLen, weight_gapO, weight_gapE, prof->profile_word, -1, maskLen);
		word = 1;
	}else {
		fprintf(stderr, "Please call the function ssw_init before ssw_align.\n");
		return 0;
	}
	r->score1 = bests[0].score;
    r->ref_end1 = bests[0].ref;
    r->read_end1 = bests[0].read;

	if (maskLen >= 15) {
		//r->score2 = bests[1].score; //mine
		//r->ref_end2 = bests[1].ref; //mine
	} else {
		//r->score2 = 0; //mine
		//r->ref_end2 = -1; //mine
	}
	free(bests);
	bests = NULL;
	if (flag == 0 || (flag == 2 && r->score1 < filters)) goto end;


	// Find the beginning position of the best alignment.
	read_reverse = seq_reverse(prof->read, r->read_end1);
	if (word == 0) {
		vP = qP_byte(read_reverse, prof->mat, r->read_end1 + 1, prof->n, prof->bias);
		bests_reverse = sw_sse2_byte(ref, 1, r->ref_end1 + 1, r->read_end1 + 1, weight_gapO, weight_gapE, vP, r->score1, prof->bias, maskLen);
	} else {
		vP = qP_word(read_reverse, prof->mat, r->read_end1 + 1, prof->n);
		bests_reverse = sw_sse2_word(ref, 1, r->ref_end1 + 1, r->read_end1 + 1, weight_gapO, weight_gapE, vP, r->score1, maskLen); //jenya
	}
	free(vP);
	vP = NULL;
	free(read_reverse);
	read_reverse = NULL;
	r->ref_begin1 = bests_reverse[0].ref;
	r->read_begin1 = r->read_end1 - bests_reverse[0].read;
    //r->read_begin1 = read_end1 - bests_reverse[0].read;
	free(bests_reverse);
	bests_reverse = NULL;
	if ((7&flag) == 0 || ((2&flag) != 0 && r->score1 < filters) || ((4&flag) != 0 && (r->ref_end1 - r->ref_begin1 > filterd || r->read_end1 - r->read_begin1 > filterd))) goto end;


	// Generate cigar.
	refLen = r->ref_end1 - r->ref_begin1 + 1;
	readLen = r->read_end1 - r->read_begin1 + 1;
	band_width = abs(refLen - readLen) + 1;
    //printf("band_width = %d\n",band_width); //TESTING
    //printf("readLen = %d\n",readLen); //TESTING
    //printf("refLen = %d\n",refLen); //TESTING
    //printf("r->score1 = %u\n",r->score1); //TESTING
    //printf("r->ref_begin1 = %d\n",r->ref_begin1); //TESTING
    //printf("r->read_begin1 = %d\n",r->read_begin1); //TESTING
    //printf("ref = %d\n",(int)*ref);
    //printf("ref-1 = %d\n",(int)*(ref-1));
    //printf("ref+r->ref_begin1 = %s\n",ref + r->ref_begin1);//TESTING
    //printf("prof->read + ->read_begin1 = %s\n",prof->read + r->read_begin1);//TESTING
    
    
	path = banded_sw(ref + r->ref_begin1, prof->read + r->read_begin1, refLen, readLen, r->score1, weight_gapO, weight_gapE, band_width, prof->mat, prof->n);


	if (path == 0) { free(r); r = NULL; } //jenya 
	else {
		r->cigar = path->seq;
		r->cigarLen = path->length;
		free(path);
		path = NULL;
	}
	
end: 

	return r;
}

void align_destroy (s_align **a) {
	if ( (*a)->cigar != 0 )
    {
		free((*a)->cigar);
        (*a)->cigar = 0;
    }
    if ( *a != 0 )
    {
        free(*a);
        *a = 0;
    }
}






/*
 * @function ssw_write() outputs the alignments reaching threshold E-value in Blast-like and/or SAM formats
 * @param char* fileout_n is the file name for the output alignments (created in paralleltraversal.cpp)
 * @param s_align *a is the alignment profile for the read (defined in ssw.h)
 * @param char* read_name is a pointer to the first character of the read name (to '@' for FASTQ and to '>' for FASTA)
 * @param char* read_seq is set to the first nucleotide in the read (if strand=0,
 *              read_seq points to the reverse-complement read)
 * @param char* read_qual points to the first character of the quality (must be set to NULL for FASTA)
 * @param char* ref_seq is in the numerical alphabet (0-4) and in the forward or 
 *              reverse-complement direction, ref_seq is set to the first nucleotide
 *              (only 1 newline at end of sequence)
 * @param double evalue is the E-value
 * @param int readlen is the read length for the current read
 * @param unsigned int bitscore is the bitscore (computed in paralleltraversal.cpp)
 * @param bool strand is set to 1 if the read is forward and to 0 if it is in reverse-complement
 * @param int blast_or_sam is set to 1 for Blast-like output and to 0 for SAM output
 */

void ssw_write (FILE** fileout,
			s_align* a,
            char* read_name,
			char* read_seq,
			char* read_qual,
			char* ref_name,	
			char* ref_seq,
			double evalue,
            int readlen,
			unsigned int bitscore,
			bool strand, // 1: forward aligned ; 0: reverse complement aligned
			int blast_or_sam) // 1: Blast like output; 0: Sam format output 
{
		char to_char[5] = {'A','C','G','T','N'};
	 
		if (blast_or_sam == 1) // Blast like output
		{
			fprintf(*fileout, "Sequence ID: ");
			char* tmp = ref_name;
			while (*tmp != '\n') fprintf(*fileout,"%c",*tmp++);
			fprintf(*fileout,"\n");

			fprintf(*fileout, "Query ID: ");
			tmp = read_name;
			while (*tmp != '\n') fprintf(*fileout,"%c",*tmp++);
			fprintf(*fileout,"\n");
			
			fprintf(*fileout, "Score: %d bits (%d)\t",bitscore,a->score1);
			fprintf(*fileout, "Expect: %.2e\t", evalue);
			if (strand) fprintf(*fileout, "strand: +\n\n");
			else fprintf(*fileout, "strand: -\n\n");
			if (a->cigar) {
				int32_t i, c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1; //mine
				while (e < a->cigarLen || left > 0) {

					int32_t count = 0;
					int32_t q = qb;
					int32_t p = pb;
					fprintf(*fileout, "Target: %8d    ", q + 1);
					for (c = e; c < a->cigarLen; ++c) {
						int32_t letter = 0xf&*(a->cigar + c);
						int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
						int32_t l = (count == 0 && left > 0) ? left: length;
						for (i = 0; i < l; ++i) {
							if (letter == 1) fprintf(*fileout, "-");
							else {
								fprintf(*fileout, "%c", to_char[(int)*(ref_seq + q)]);
								++q;
							}
							++ count;
							if (count == 60) goto step2;
						}
					}
	step2:
					fprintf(*fileout, "    %d\n                    ", q);
					q = qb;
					count = 0;
					for (c = e; c < a->cigarLen; ++c) {
						int32_t letter = 0xf&*(a->cigar + c);
						int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
						int32_t l = (count == 0 && left > 0) ? left: length;
						for (i = 0; i < l; ++i){ 
							if (letter == 0) {
								if ((char)to_char[(int)*(ref_seq + q)] == (char)to_char[(int)*(read_seq + p)]) fprintf(*fileout, "|");
								else fprintf(*fileout, "*");
								++q;
								++p;
							} else {
								fprintf(*fileout, " ");
								if (letter == 1) ++p;
								else ++q;
							}
							++ count;
							if (count == 60) {
								qb = q;
								goto step3;
							}
						}
					}
	step3:
					p = pb;
					fprintf(*fileout, "\nQuery:  %8d    ", p + 1);
					count = 0;
					for (c = e; c < a->cigarLen; ++c) {
						int32_t letter = 0xf&*(a->cigar + c);
						int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
						int32_t l = (count == 0 && left > 0) ? left: length;
						for (i = 0; i < l; ++i) { 
							if (letter == 2) fprintf(*fileout, "-");
							else {							
								fprintf(*fileout, "%c", (char)to_char[(int)*(read_seq + p)]);
								++p;
							}
							++ count;
							if (count == 60) {
								pb = p;
								left = l - i - 1;
								e = (left == 0) ? (c + 1) : c;
								goto end;
							}
						}
					}
					e = c;
					left = 0;
	end:
					fprintf(*fileout, "    %d\n\n", p);
				}
			}
		}
		/* SAM output */
		else
		{
			/// (1) QNAME 
			while ((*read_name != ' ') && (*read_name != '\n') && (*read_name != '\t')) fprintf(*fileout,"%c",*read_name++);
            /* set the read name ptr to end of read name */
            while (*read_name != '\n') read_name++;

			if (a->score1 == 0) fprintf(*fileout, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
			else {
				int32_t c;
                //int32_t l = a->read_end1 - a->read_begin1 + 1;
                int32_t qb = a->ref_begin1;
                int32_t pb = a->read_begin1;
                int32_t p;
				//uint32_t mapq = -4.343 * log(1 - (double)abs(a->score1 - a->score2)/(double)a->score1); //jenya
				//mapq = (uint32_t) (mapq + 4.99); //jenya
				//mapq = mapq < 254 ? mapq : 254; //jenya

				/* (2) FLAG */
				if (!strand) fprintf(*fileout, "\t16\t");
				else fprintf(*fileout, "\t0\t");

				/* (3) RNAME */
				while ((*ref_name != ' ') && (*ref_name != '\n') && (*ref_name != '\t')) fprintf(*fileout,"%c",*ref_name++);

				/* (4) POS */
				fprintf(*fileout, "\t%d", a->ref_begin1+1);

				/* (5) MAPQ */
				fprintf(*fileout,"\t%d\t", 255);

				/* (6) CIGAR */
                if ( a->read_begin1 != 0 ) //output the masked region at beginning of alignment
                    fprintf(*fileout,"%dS",a->read_begin1);
                
                int32_t len_align = 0; //TESTING
                
				for (c = 0; c < a->cigarLen; ++c) {
					int32_t letter = 0xf&*(a->cigar + c);
					int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
                    //if ( letter != 2 ) len_align+=length; //TESTING
					fprintf(*fileout, "%d", length);
					if (letter == 0) { fprintf(*fileout, "M"); len_align+=length;}
					else if (letter == 1) { fprintf(*fileout, "I"); len_align+=length;}
					else fprintf(*fileout, "D");
				}
                
                //int32_t end_mask = readlen-((a->read_end1)-(a->read_begin1)+1)-(a->read_begin1);
                int32_t end_mask = readlen-(a->read_begin1)-len_align;
                if ( end_mask > 0 ) //output the masked region at end of alignment
                    fprintf(*fileout,"%dS",end_mask);
                
                //printf("a->read_begin1 = %d\n",a->read_begin1); //TESTING
                //printf("a->read_end1 = %d\n",a->read_end1); //TESTING
                //printf("readlen = %d\n", readlen); //TESTING
                //printf("cigarlen = %d\n",len_align); //TESTING
                //printf("end_mask = %d\n\n",end_mask); //TESTING
                //printf("read_end1 - read_begin1 = %d\n",((a->read_end1)-(a->read_begin1)+1)); //TESTING
                //printf("cigar length = %d\n\n",len_align); //TESTING
                
				/* (7) RNEXT, (8) PNEXT, (9) TLEN */
				fprintf(*fileout, "\t*\t0\t0\t");

				/* (10) SEQ */
                char* ptr_read_seq = read_seq;
                while (*ptr_read_seq != '\n') fprintf(*fileout,"%c",(char)to_char[(int)*ptr_read_seq++]);
                                
				/* (11) QUAL */
                fprintf(*fileout,"\t");

                /* reverse-complement strand */
				if (read_qual && !strand)
                {
                   while (*read_qual != '\n') fprintf(*fileout,"%c",*read_qual--);
                /* forward strand */
				}else if (read_qual){
                    while ((*read_qual != '\n') && (*read_qual!='\0')) fprintf(*fileout,"%c",*read_qual++);
                /* FASTA read */
				} else fprintf(*fileout, "*");


				/* (12) OPTIONAL FIELD: SW alignment score generated by aligner */
				fprintf(*fileout, "\tAS:i:%d", a->score1);

				/* (13) OPTIONAL FIELD: edit distance to the reference */
                uint32_t diff = 0;
                for (c = 0; c < a->cigarLen; ++c) {
                    int32_t letter = 0xf&*(a->cigar + c);
                    int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
                    if (letter == 0) {
                        for (p = 0; p < length; ++p){
                            if ( (char)to_char[(int)*(ref_seq + qb)] != (char)*(read_seq + pb) ) ++diff;
                            ++qb;
                            ++pb;
                    }
                    } else if (letter == 1) {
                        pb += length;
                        diff += length;
                    } else {
                        qb += length;
                        diff += length;
                    }
                }

				fprintf(*fileout,"\tNM:i:%d\n", diff);

			}
		}//~sam output
}



/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

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

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef BWA_BWT_H
#define BWA_BWT_H

#include <stdint.h>
#include <stddef.h>

// requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code assume OCC_INTERVAL=0x80
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef unsigned char ubyte_t;
#endif

typedef struct {
    uint64_t *bwt0;
    uint64_t *bwt1;
    size_t size;
} bwt_vector;

typedef struct {
    uint64_t *occurrences;
    size_t rows;
    size_t cols;
} bwt_occ_matrix;

typedef uint64_t bwtint_t; //unsigned long long 2^64

//gives name to user defined data type --> bwt_t
typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	//irst occurrence of A..., C..., G..., T.... and $
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint64_t *bwt_new;

	// keep this constant
	// occurrence array, separated to two parts
	uint32_t cnt_table[256];

    bwt_vector subseq;

	// suffix array
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
    uint8_t parent;               /* parent or daughter */

    //char before $
    int c;

    //new data structure
    bwtint_t bwt0_size;

} bwt_t;

/**
 * x[0] - forward index location;
 * x[1] - reverse complement index location;
 * x[2] - occurrence number (number of substrings);
 * info>>32 - beg;
 * (uint32_t) info - end;
 */

//gives name to user defined data type --> bwtintv_t
typedef struct {
	bwtint_t x[3], info;
} bwtintv_t;

typedef struct {
    size_t n, m;
    bwtintv_t *a;
} bwtintv_v;

/* For general OCC_INTERVAL, the following is correct:
#define bwt_bwt(b, k) ((b)->bwt[(k)/OCC_INTERVAL * (OCC_INTERVAL/(sizeof(uint32_t)*8/2) + sizeof(bwtint_t)/4*4) + sizeof(bwtint_t)/4*4 + (k)%OCC_INTERVAL/16])
#define bwt_occ_intv(b, k) ((b)->bwt + (k)/OCC_INTERVAL * (OCC_INTERVAL/(sizeof(uint32_t)*8/2) + sizeof(bwtint_t)/4*4)
*/

// The following two lines are ONLY correct when OCC_INTERVAL==0x80
#define bwt_bwt(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))

/* retrieve a character from the $-removed BWT string. Note that
 * bwt_t::bwt is not exactly the BWT string and therefore this macro is
 * called bwt_B0 instead of bwt_B */
#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

/* #define bwt_set_intv(bwt, c, ik) ((ik).x[0] = (bwt)->L2[(int)(c)]+1, (ik).x[2] = (bwt)->L2[(int)(c)+1]-(bwt)->L2[(int)(c)], (ik).x[1] = (bwt)->L2[3-(c)]+1, (ik).info = 0) */
// interval for the first base
#define bwt_set_intv(bwt, bwtc, c, ik) ((ik).x[0] = (bwt)->L2[(int)(c)]+1, (ik).x[2] = (bwt)->L2[(int)(c)+1]-(bwt)->L2[(int)(c)], (ik).x[1] = (bwtc)->L2[3-(c)]+1, (ik).info = 0)

#ifdef __cplusplus
extern "C" {
#endif

	void bwt_dump_bwt(const char *fn, const bwt_t *bwt);
	void bwt_dump_sa(const char *fn, const bwt_t *bwt);

	bwt_t *bwt_restore_bwt(const char *fn);
  void bwt_restore_bwt2(const char *fn, bwt_t *bwt);
	void bwt_restore_sa(const char *fn, bwt_t *bwt);

	void bwt_destroy(bwt_t *bwt);

	void bwt_bwtgen(const char *fn_pac, const char *fn_bwt); // from BWT-SW
	void bwt_bwtgen2(const char *fn_pac, const char *fn_bwt, int block_size); // from BWT-SW
	void bwt_cal_sa(bwt_t *bwt, int intv);

	void bwt_bwtupdate_core(bwt_t *bwt, int index);

	bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, int c);
	void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4]);
	bwtint_t bwt_sa(bwt_t *bwt, bwtint_t k);

	// more efficient version of bwt_occ/bwt_occ4 for retrieving two close Occ values
	void bwt_gen_cnt_table(bwt_t *bwt);
	void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
	void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4]);

	int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end);
	int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0);

	/**
	 * Extend bi-SA-interval _ik_
	 */
	void bwt_extend(const bwt_t *bwt, bwtintv_t *ik, bwtintv_t *ok, int is_back, int c);

	/**
	 * Given a query _q_, collect potential SMEMs covering position _x_ and store them in _mem_.
	 * Return the end of the longest exact match starting from _x_.
	 */
	int bwt_smem1(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]);
	int bwt_smem1a(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]);

	int bwt_seed_strategy1(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem);

#ifdef __cplusplus
}
#endif

#endif

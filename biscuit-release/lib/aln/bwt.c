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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include "utils.h"
#include "bwt.h"
#include "kvec.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

void bwt_gen_cnt_table(bwt_t *bwt) {
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k) // compute inverse CSA
{
	//fprintf(stderr, "k: %llu\n", k);
	bwtint_t x = k - (k > bwt->primary);
	//fprintf(stderr, "x: %llu\n", x);
	x = bwt_B0(bwt, x);
	//fprintf(stderr, "x: %llu\n", x);
	//fprintf(stderr, "bwt->L2[x]: %llu\n", bwt->L2[x]);
	//fprintf(stderr, "bwt_occ(bwt, k, x): %llu\n", bwt_occ(bwt, k, x));
	x = bwt->L2[x] + bwt_occ(bwt, k, x);
	//fprintf(stderr, "x: %llu\n", x);
	return k == bwt->primary? 0 : x;
}

// bwt->bwt and bwt->occ must be precalculated
void bwt_cal_sa(bwt_t *bwt, int intv)
{
	bwtint_t isa, sa, i; // S(isa) = sa
	int intv_round = intv;

	kv_roundup32(intv_round);
	xassert(intv_round == intv, "SA sample interval is not a power of 2.");
	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len/2 + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	// calculate SA value
	isa = 0; sa = bwt->seq_len/2;
	for (i = 0; i < bwt->seq_len/2; ++i) {

		if (isa % intv == 0) {
		    bwt->sa[isa/intv] = sa;
		    fprintf(stderr, "sa: %llu isa: %llu\n", sa, isa);
		}
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;
	bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t sa = 0, mask = bwt->sa_intv - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv];
}

static inline int __occ_aux(uint64_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
	bwtint_t n;
	uint32_t *p, *end;

	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (bwtint_t)(-1)) return 0;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// retrieve Occ at k/OCC_INTERVAL
	n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
	p += sizeof(bwtint_t); // jump to the start of the first BWT cell

	// calculate Occ up to the last k/32
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

	// calculate Occ
	n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31; // corrected for the masked bits

	return n;
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
{
	bwtint_t _k, _l;
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		*ok = bwt_occ(bwt, k, c);
		*ol = bwt_occ(bwt, l, c);
	} else {
		bwtint_t m, n, i, j;
		uint32_t *p;
		if (k >= bwt->primary) --k;
		if (l >= bwt->primary) --l;
		n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
		p += sizeof(bwtint_t);
		// calculate *ok
		j = k >> 5 << 5;
		for (i = k/OCC_INTERVAL*OCC_INTERVAL; i < j; i += 32, p += 2)
			n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
		m = n;
		n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
		if (c == 0) n -= ~k&31; // corrected for the masked bits
		*ok = n;
		// calculate *ol
		j = l >> 5 << 5;
		for (; i < j; i += 32, p += 2)
			m += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
		m += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
		if (c == 0) m -= ~l&31; // corrected for the masked bits
		*ol = m;
	}
}

#define __occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

// and here is how to find the occurrence of k not divisable by 128.
void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	bwtint_t x;
	uint32_t *p, tmp, *end;
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	k -= (k >= bwt->primary); // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	memcpy(cnt, p, 4 * sizeof(bwtint_t));
	p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
	end = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4)); // this is the end point of the following loop
	for (x = 0; p < end; ++p) x += __occ_aux4(bwt, *p);
	tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
	x += __occ_aux4(bwt, tmp) - (~k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}
//for inverse, when seq0 and seq1 ^= ULLONG_MAX
//after adjusted for k
void __builtin_popcountll_inv_occ(uint64_t seq0, uint64_t seq1, bwtint_t counts[3]) {

    //fprintf(stderr,"bwt0: %llu bwt1: %llu\n", seq0, seq1);
    //fprintf(stderr,"counts[0]: %llu counts[1]: %llu counts[2]: %llu\n", counts[0], counts[1], counts[2]);

    uint64_t partial_G = __builtin_popcountll(seq0 & seq1);
    //G
    counts[2] += partial_G;
    //A
    counts[0] += __builtin_popcountll(seq1) - partial_G;
    //T
    counts[1] += __builtin_popcountll(seq0) - partial_G;

    //fprintf(stderr,"counts[0]: %llu counts[1]: %llu counts[2]: %llu\n", counts[0], counts[1], counts[2]);
}

//k ranking will start at 1; can be changed later to start at 0
void bwt_occ4_new_index(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[3])
{
//fix to take into value of k

    int fix = 0;

    //total values
    if ((double) k / 128 > 1) {
        cnt[0] = bwt->bwt_occ_matrix0.occurrences[(k/128 - 1) * 3 + 0];
        cnt[1] = bwt->bwt_occ_matrix0.occurrences[(k/128 - 1) * 3 + 1];
        cnt[2] = bwt->bwt_occ_matrix0.occurrences[(k/128 - 1) * 3 + 2];
    } else {
        cnt[0] = cnt[1] = cnt[2] = 0;
    }

    //we might need this again
    if ((double) (k % 128) / 64 > 1)
        fix = 1;

    //fprintf(stderr, "At last checkpoint ->\nA: %llu T: %llu G: %llu\n", cnt[0], cnt[1], cnt[2]);

    bwt_t *bwt_k;

    bwt_k = (bwt_t*)calloc(1, sizeof(bwt_t));

    uint8_t l = ceil((double)(k % 128)/64);

    // intialize space of vector
    bwt_k->bwt_all.bwt0 = (uint64_t*)calloc(l, 8);
    bwt_k->bwt_all.bwt1 = (uint64_t*)calloc(l, 8);

    for (int i = 0; i < l; i++) {

        bwt_k->bwt_all.bwt0[i] = bwt->bwt_all.bwt0[(int)floor(k/64) + i - fix];
        bwt_k->bwt_all.bwt1[i] = bwt->bwt_all.bwt1[(int)floor(k/64) + i - fix];

        //fprintf(stderr, "BWT0: %llu BWT1: %llu\n", bwt_k->bwt_all.bwt0[0], bwt_k->bwt_all.bwt1[i]);

        bwt_k->bwt_all.bwt0[i] ^= ULLONG_MAX;
        bwt_k->bwt_all.bwt1[i] ^= ULLONG_MAX;

        //this k adjustment needs to be done before calling __builtin_popcountll_inv_occ
        if (i == l - 1) {
            bwt_k->bwt_all.bwt0[i] = bwt_k->bwt_all.bwt0[i] >> (64 - k % 64);
            bwt_k->bwt_all.bwt1[i] = bwt_k->bwt_all.bwt1[i] >> (64 - k % 64);
        }

        __builtin_popcountll_inv_occ(bwt_k->bwt_all.bwt0[i], bwt_k->bwt_all.bwt1[i], cnt);

    }

    fprintf(stderr, "Up to position %llu -->\nA: %llu T: %llu G: %llu\n", k, cnt[0], cnt[1], cnt[2]);

    //free(bwt_k);
    //free(bwt_k->bwt_all.bwt0);
    //free(bwt_k->bwt_all.bwt1);

}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	//bwtint_t _k, _l;
	//_k = k - (k >= bwt->primary);
	//_l = l - (l >= bwt->primary);

    bwt_occ4_new_index(bwt, k, cntk);
    bwt_occ4_new_index(bwt, l, cntl);
	/*
	if (_l>>OCC_INTV_SHIFT != _k>>OCC_INTV_SHIFT || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4_new_index(bwt, k, cntk);
		bwt_occ4_new_index(bwt, l, cntl);
		fprintf(stderr, "TEST\n");
	} else {
		bwtint_t x, y;
		uint32_t *p, tmp, *endk, *endl;
		k -= (k >= bwt->primary); // because $ is not in bwt
		l -= (l >= bwt->primary);
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
		// prepare cntk[]
		endk = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4));
		endl = p + ((l>>4) - ((l&~OCC_INTV_MASK)>>4));

		//is it okay to replace this with the new index? not sure what __occ_aux4
		//or how the cnt_table array is formatted
		for (x = 0; p < endk; ++p)
		    x += __occ_aux4(bwt, *p);
		y = x;
		tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
		x += __occ_aux4(bwt, tmp) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		for (; p < endl; ++p) y += __occ_aux4(bwt, *p);
		tmp = *p & ~((1U<<((~l&15)<<1)) - 1);
		y += __occ_aux4(bwt, tmp) - (~l&15);
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
	*/
}

int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t k, l, ok, ol;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) break; // no match
	}
	if (k > l) return 0; // no match
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
	int i;
	bwtint_t k, l, ok, ol;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}

/*********************
 * Bidirectional BWT *
 *********************/

void bwt_extend(const bwt_t *bwt, bwtintv_t *ik, bwtintv_t ok[3], int is_back, uint64_t *sub, int size) //add argument for the sequence
{
    bwtint_t total[3];
    bwt_occ4_new_index(bwt, bwt->seq_len/2, total);

    uint64_t total_A = total[0];
    uint64_t total_G = total[2];


    uint64_t *sub_ = (uint64_t*)calloc(2, 8);
    bwtint_t cntk_old[3];
    bwtint_t cntk_new[3];
    int first;
    cntk_old[0] = cntk_new[0] = 0;    cntk_old[1] = cntk_new[1] = 0;    cntk_old[2] = cntk_new[2] = 0;

    //pre adjustment for __builtin_popcountll_inv_occ
    sub[0] ^= ULLONG_MAX;
    sub[1] ^= ULLONG_MAX;
    sub_[0] = sub[0] >> (64 - size);
    sub_[1] = sub[1] >> (64 - size);

    __builtin_popcountll_inv_occ(sub_[0], sub_[1], cntk_old);
    fprintf(stderr, "Composition at pos %d -->\nA: %llu T: %llu G: %llu\n", size - 1,
        cntk_old[0], cntk_old[1], cntk_old[2]);
    sub_[0] = sub[0] >> (64 - size + 1);
    sub_[1] = sub[1] >> (64 - size + 1);
    __builtin_popcountll_inv_occ(sub_[0], sub_[1], cntk_new);
    fprintf(stderr, "Nucleotide at pos %d -->\nA: %llu T: %llu G: %llu\n\n",
        size - 1, cntk_old[0] - cntk_new[0], cntk_old[1] - cntk_new[1], cntk_old[2] - cntk_new[2]);

    /*     0 = A
    *     2 = G
    *     3 = T
    */

    if (cntk_old[0] - cntk_new[0] == 1)
        first = 0;
    else if (cntk_old[1] - cntk_new[1] == 1)
        first = 3;
    else if (cntk_old[2] - cntk_new[2])
        first = 2;

    bwt_set_intv(bwt, bwt, first, *ik);

	bwtint_t tk[3], tl[3];
	int i;

    uint64_t k = ik->x[!is_back] - 1;
    uint64_t l = ik->x[!is_back] - 1 + ik->x[2];

    //set the first value
    for (int j = 1; j < size + 1; j++) {

        fprintf(stderr, "k: %d l: %d\n", k, l);

        if (j == size) {
            break;
        }

        bwt_2occ4(bwt, k, l, tk, tl);

        //A
        fprintf(stderr, "tk[0]: %llu tl[0]: %llu\n", tk[0], tl[0]);
        //T
        fprintf(stderr, "tk[1]: %llu tl[1]: %llu\n", tk[1], tl[1]);
        //G
        fprintf(stderr, "tk[2]: %llu tl[2]: %llu\n", tk[2], tl[2]);

        cntk_old[0] = cntk_new[0];  cntk_old[1] = cntk_new[1];  cntk_old[2] = cntk_new[2];

        fprintf(stderr, "Composition at pos %d -->\nA: %llu T: %llu G: %llu\n", size - j,
            cntk_old[0], cntk_old[1], cntk_old[2]);

        cntk_new[0] = cntk_new[1] = cntk_new[2] = 0;

        sub_[0] = sub[0] >> (64 - size + 1 + j);
        sub_[1] = sub[1] >> (64 - size + 1 + j);
        if (j != size - 1) {
                __builtin_popcountll_inv_occ(sub_[0], sub_[1], cntk_new);
        } else {
            //cntk_new[0], cntk_new[1] = cntk_new[2] = 0;
            //this needs to be the output of the function
            //work on this next
            fprintf(stderr, "k: %d l: %d\n", k, l);
            //fprintf(stderr);
        }

        fprintf(stderr, "Composition at pos %d -->\nA: %llu T: %llu G: %llu\n", size - j - 1,
            cntk_new[0], cntk_new[1], cntk_new[2]);

        fprintf(stderr, "Nucleotide at pos %d -->\nA: %llu T: %llu G: %llu\n\n",
            size - 1 - j, cntk_old[0] - cntk_new[0], cntk_old[1] - cntk_new[1], cntk_old[2] - cntk_new[2]);

        //cnt 0 = A     1 = T       2 = G

        if (cntk_old[0] - cntk_new[0] == 1) {
            k = tk[0];
            l = tl[0];
            fprintf(stderr, "A\n");
        } else if (cntk_old[1] - cntk_new[1] == 1) {
            k = total_A + total_G + tk[1];
            l = total_A + total_G + tl[1];
        } else if (cntk_old[2] - cntk_new[2] == 1) {
            k = total_A + tk[2];
            l = total_A + tl[2];
        } else if (cntk_old[0] == cntk_new[0] && cntk_old[1] == cntk_new[1] &&
            cntk_old[2] == cntk_new[2]) {
            break;
        }
        //tk[0] = tk[1] = tk[2] = 0;
        //tl[0] = tl[1] = tk[2] = 0;
    }

exit(0);
    //change so that
	for (i = 0; i != 3; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}



	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
	ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
	ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];

	fprintf(stderr, "ok[3].x[is_back]: %llu\n", ok[3].x[is_back]);
	fprintf(stderr, "ok[2].x[is_back]: %llu\n", ok[2].x[is_back]);
	fprintf(stderr, "ok[1].x[is_back]: %llu\n", ok[1].x[is_back]);
	fprintf(stderr, "ok[0].x[is_back]: %llu\n", ok[0].x[is_back]);

}

static void bwt_reverse_intvs(bwtintv_v *p) {
  if (p->n > 1) {
    int j;
    for (j = 0; (unsigned) j < p->n>>1; ++j) {
      bwtintv_t tmp = p->a[p->n - 1 - j];
      p->a[p->n - 1 - j] = p->a[j];
      p->a[j] = tmp;
    }
  }
}
// NOTE: $max_intv is not currently used in BWA-MEM
int bwt_smem1a(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]) {

  int i, j, c, ret;
  bwtintv_t ik, ok[4];
  bwtintv_v a[2], *prev, *curr, *swap;

  mem->n = 0;
  if (q[x] > 3) return x + 1;
  if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
  kv_init(a[0]); kv_init(a[1]);
  prev = tmpvec && tmpvec[0]? tmpvec[0] : &a[0]; // use the temporary vector if provided
  curr = tmpvec && tmpvec[1]? tmpvec[1] : &a[1];
  bwt_set_intv(bwt, bwtc, q[x], ik); // the initial interval of a single base
  ik.info = x + 1;

  for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
    if (ik.x[2] < max_intv) { // an interval small enough
      kv_push(bwtintv_t, *curr, ik);
      break;
    } else if (q[i] < 4) { // an A/C/G/T base
      c = 3 - q[i]; // complement of q[i]
      //bwt_extend(bwtc, &ik, ok, 0);
      if (ok[c].x[2] != ik.x[2]) { // change of the interval size
        kv_push(bwtintv_t, *curr, ik);
        if (ok[c].x[2] < (unsigned) min_intv) break; // the interval size is too small to be extended further
      }
      ik = ok[c]; ik.info = i + 1;
    } else { // an ambiguous base
      kv_push(bwtintv_t, *curr, ik);
      break; // always terminate extension at an ambiguous base; in this case, i<len always stands
    }
  }
  if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
  bwt_reverse_intvs(curr); // s.t. smaller intervals (i.e. longer matches) visited first
  ret = curr->a[0].info; // this will be the returned value
  swap = curr; curr = prev; prev = swap;

  for (i = x - 1; i >= -1; --i) { // backward search for MEMs
    c = i < 0? -1 : q[i] < 4? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
    for (j = 0, curr->n = 0; (unsigned) j < prev->n; ++j) {
      bwtintv_t *p = &prev->a[j];
      //if (c >= 0 && ik.x[2] >= max_intv) bwt_extend(bwt, p, ok, 1);
      if (c < 0 || ik.x[2] < (unsigned) max_intv || ok[c].x[2] < (unsigned) min_intv) { // keep the hit if reaching the beginning or an ambiguous base or the intv is small enough
        if (curr->n == 0) { // test curr->n>0 to make sure there are no longer matches
          if (mem->n == 0 || (unsigned) i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
            ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
            kv_push(bwtintv_t, *mem, ik);
          }
        } // otherwise the match is contained in another longer match
      } else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]) {
        ok[c].info = p->info;
        kv_push(bwtintv_t, *curr, ok[c]);
      }
    }
    if (curr->n == 0) break;
    swap = curr; curr = prev; prev = swap;
  }
  bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

  if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
  if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
  return ret;
}

int bwt_smem1(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]) {
  return bwt_smem1a(bwt, bwtc, len, q, x, min_intv, 0, mem, tmpvec);
}

int bwt_seed_strategy1(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem) {
  int i, c;
  bwtintv_t ik, ok[4];

  memset(mem, 0, sizeof(bwtintv_t));
  if (q[x] > 3) return x + 1;
  bwt_set_intv(bwt, bwtc, q[x], ik); // the initial interval of a single base
  for (i = x + 1; i < len; ++i) { // forward search
    if (q[i] < 4) { // an A/C/G/T base
      c = 3 - q[i]; // complement of q[i]
      //bwt_extend(bwtc, &ik, ok, 0);
      if (ok[c].x[2] < (unsigned) max_intv && i - x >= min_len) {
        *mem = ok[c];
        mem->info = (uint64_t)x<<32 | (i + 1);
        return i + 1;
      }
      ik = ok[c];
    } else return i + 1;
  }
  return len;
}

/*************************
 * Read/write BWT and SA *
 *************************/

void bwt_dump_bwt(const char *fn, const bwt_t *bwt) {
  FILE *fp;
  fp = xopen(fn, "wb");
  err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
  err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
  err_fwrite(bwt->bwt, 4, bwt->bwt_size, fp);
  err_fflush(fp);
  err_fclose(fp);
}

void bwt_dump_sa(const char *fn, const bwt_t *bwt) {
  FILE *fp;
  fp = xopen(fn, "wb");
  err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
  err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
  err_fwrite(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
  err_fwrite(&bwt->seq_len, sizeof(bwtint_t), 1, fp);
  err_fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
  err_fflush(fp);
  err_fclose(fp);
}

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a) { // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
  const int bufsize = 0x1000000; // 16M block
  bwtint_t offset = 0;
  while (size) {
    int x = (unsigned) bufsize < size ? bufsize : (int) size;
    if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
    size -= x; offset += x;
  }
  return offset;
}

void bwt_restore_sa(const char *fn, bwt_t *bwt) {
  char skipped[256];
  FILE *fp;
  bwtint_t primary;

  fp = xopen(fn, "rb");
  err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
  xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
  err_fread_noeof(skipped, sizeof(bwtint_t), 4, fp); // skip
  err_fread_noeof(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
  err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
  xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

  bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
  bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
  bwt->sa[0] = -1;

  fread_fix(fp, sizeof(bwtint_t) * (bwt->n_sa - 1), bwt->sa + 1);
  err_fclose(fp);
}

bwt_t *bwt_restore_bwt(const char *fn) {
  bwt_t *bwt;
  FILE *fp;

  bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
  fp = xopen(fn, "rb");
  err_fseek(fp, 0, SEEK_END);
  bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
  bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
  err_fseek(fp, 0, SEEK_SET);
  err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
  err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
  fread_fix(fp, bwt->bwt_size<<2, bwt->bwt);
  bwt->seq_len = bwt->L2[4];
  err_fclose(fp);
  bwt_gen_cnt_table(bwt);

  return bwt;
}

void bwt_restore_bwt2(const char *fn, bwt_t *bwt) {
  FILE *fp;
  memset(bwt, 0, sizeof(bwt_t));
  /* bwt = (bwt_t*)calloc(1, sizeof(bwt_t)); */
  fp = xopen(fn, "rb");
  err_fseek(fp, 0, SEEK_END);
  bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
  bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
  err_fseek(fp, 0, SEEK_SET);
  err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
  err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
  fread_fix(fp, bwt->bwt_size<<2, bwt->bwt);
  bwt->seq_len = bwt->L2[4];
  err_fclose(fp);
  bwt_gen_cnt_table(bwt);
}

void bwt_destroy(bwt_t *bwt) {
  if (bwt == 0) return;
  free(bwt->sa); free(bwt->bwt);
  free(bwt);
}



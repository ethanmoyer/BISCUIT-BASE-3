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

// This function is called to get the counts for the first k nucleotides. Here k is accepted as positioning starting at
// 1.
bwtint_t builtin_popcountll(uint64_t seq0, uint64_t seq1, int c, bwtint_t k) {
    uint64_t partial_G = 0;
    uint64_t top = seq0 >> (64 - k);
    uint64_t bottom = seq1 >> (64 - k);
    switch(c) {
        case 0:
            return __builtin_popcountll(top & ~bottom); //A
        case 1:
            return __builtin_popcountll(top & bottom); //C
        case 2:
            partial_G = __builtin_popcountll(~top & ~bottom);
            return k % 64 != 0 ? partial_G - 64 + k % 64 : partial_G; //G
        case 3:
            return __builtin_popcountll(~top & bottom); //T
    }
    return 0;
}
// This is the bwt_occ equivalent for gathering counts when a position isn't divisible by 128.

bwtint_t bwt_occ_new_index(const bwt_t *bwt, bwtint_t k, int c) {
    // bwt starts indexing at 0, so calculations involving k are handled accordingly.
    k -= (k >= bwt->primary);

    bwtint_t count = 0;

    if ((k + 1) / 128)
        count = bwt->bwt_new[((k + 1) / 128 - 1) * 8 + c];

    if ((k + 1) == bwt->seq_len && c == 2 && count == 0)
        return 0;

    if ((k + 1) % 128 == 0) return count;

    if ((k + 1) % 128 < 64) {
        count += builtin_popcountll(bwt->bwt_new[k/128 * 8 + 4], bwt->bwt_new[k/128 * 8 + 6], c, ((k + 1) % 128));
    } else {
        count += builtin_popcountll(bwt->bwt_new[k/128 * 8 + 4], bwt->bwt_new[k/128 * 8 + 6], c, 64);
        if ((k + 1) % 64 == 0) return count;
        count += builtin_popcountll(bwt->bwt_new[k/128 * 8 + 5], bwt->bwt_new[k/128 * 8 + 7], c, (k % 64) + 1);
    }
    return count;
}

// Retrieves the base at a specific position with the least amount of operations as possible (or the least
// to my knowledge). There may be room for improvement.
// Seg fault is here
int nucAtBWTinv(bwt_t *bwt, bwtint_t k) {
    bwtint_t top = (bwt->bwt_new[k/128 * 8 + ((k % 128) < 64 ? 4 : 5)] >> (63 - k % 64)) & 1;
    bwtint_t bottom = (bwt->bwt_new[k/128 * 8 + ((k % 128) < 64 ? 6 : 7)] >> (63 - k % 64)) & 1;
    bwtint_t p = top | bottom;
    if (p == 0)
        return 2; //return G
    else if (~top & bottom)
        return 3; //return T
    else if (top & ~bottom)
        return 0; //return A
    return p; //return C
}

// compute inverse CSA
static inline bwtint_t bwt_invPsi(bwt_t *bwt, bwtint_t k) {
    bwtint_t x = k - (k > bwt->primary);

    x = nucAtBWTinv(bwt, x);
    x = bwt->L2[x] + bwt_occ_new_index(bwt, k, x);
    return k == bwt->primary ? 0 : x;
}

void bwt_cal_sa(bwt_t *bwt, int intv) {
    bwtint_t isa, sa, i; // S(isa) = sa
    int intv_round = intv;
    kv_roundup32(intv_round);
    xassert(intv_round == intv, "SA sample interval is not a power of 2.");
    xassert(bwt->bwt_new, "bwt_t::bwt is not initialized.");
    if (bwt->sa) free(bwt->sa);
    bwt->sa_intv = intv;
    bwt->n_sa = (bwt->seq_len + intv) / intv;
    bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
    // calculate SA value
    isa = 0; sa = bwt->seq_len;
    for (i = 0; i < bwt->seq_len; ++i) {
        if (isa % intv == 0) {
            bwt->sa[isa/intv] = sa;
        }
        --sa;
        isa = bwt_invPsi(bwt, isa);
    }
    if (isa % intv == 0) bwt->sa[isa/intv] = sa;
    //bwt->sa[0] = 0; //might change back later
    //bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

bwtint_t bwt_sa(bwt_t *bwt, bwtint_t k) {
    bwtint_t sa = 0, mask = bwt->sa_intv - 1;
    while (k % mask) { // != 0
        ++sa;
        k = bwt_invPsi(bwt, k);
    }
    /* without setting bwt->sa[0] = -1, the following line should be
       changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
    return sa + bwt->sa[k/bwt->sa_intv];
}

void bwt_gen_cnt_table(bwt_t *bwt) {
    int i, j;
    for (i = 0; i != 256; ++i) {
        uint32_t x = 0;
        for (j = 0; j != 4; ++j)
            x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
        bwt->cnt_table[i] = x;
    }
}

static inline int occ_aux(uint64_t y, int c)
{
    // reduce nucleotide counting to bits counting
    y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
    // count the number of 1s in y
    y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
    return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, int c)
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
    for (; p < end; p += 2) n += occ_aux((uint64_t)p[0]<<32 | p[1], c);

    // calculate Occ
    n += occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
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
            n += occ_aux((uint64_t)p[0]<<32 | p[1], c);
        m = n;
        n += occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
        if (c == 0) n -= ~k&31; // corrected for the masked bits
        *ok = n;
        // calculate *ol
        j = l >> 5 << 5;
        for (; i < j; i += 32, p += 2)
            m += occ_aux((uint64_t)p[0]<<32 | p[1], c);
        m += occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
        if (c == 0) m -= ~l&31; // corrected for the masked bits
        *ol = m;
    }
}

#define occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

//return four
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
	for (x = 0; p < end; ++p) x += occ_aux4(bwt, *p);
	tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
	x += occ_aux4(bwt, tmp) - (~k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
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

void bwt_extend(const bwt_t *bwt, bwtintv_t *ik, bwtintv_t *ok, int is_back, int c) {
    bwtint_t tk[4], tl[4];
    //bwt_occ_new_index_v2(bwt, ik->x[!is_back] - 1, ik->x[!is_back] + ik->x[2] - 1, tk, tl, c);
    for (int i = 3; i != c - 1; i--) {
        tk[i] = bwt_occ_new_index(bwt, ik->x[!is_back] - 1, i);
        tl[i] = bwt_occ_new_index(bwt, ik->x[!is_back] + ik->x[2] - 1, i);
        ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
        ok[i].x[2] = tl[i] - tk[i];
        if (i == 3)
            ok[3].x[is_back] = ik->x[is_back] + (
                    ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
        ok[i - 1].x[is_back] = ok[i].x[is_back] + ok[i].x[2];
    }
    // Implemented it all in one for loop... will need to test.
    // if i starts at c in bwt_occ_new_index_v2, then go from 3 to c in the stuff below
    //ok[3].x[is_back] = ik->x[is_back] + (
    //      ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
    //for (int i = 3; i > c; --i) {
    //}

    //ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
    //ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
    //ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
}

int bwt_smem1a (const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]) {

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
        if (ik.x[2] < max_intv) { // an interval small enough, currently won't come here, max_intv is currently 0
            kv_push(bwtintv_t, *curr, ik);
            break;
        } else if (q[i] < 4) { // an A/C/G/T base
            c = 3 - q[i];
            bwt_extend(bwtc, &ik, ok, 0, c);
            if (ok[c].x[2] != ik.x[2]) { // change of the interval size
                kv_push(bwtintv_t, *curr, ik);
                if (ok[c].x[2] < (unsigned) min_intv) { // no more matches
                    break; // the interval size is too small to be extended further
                }
            }
            ik = ok[c]; ik.info = i + 1;
        } else { // an ambiguous base
            kv_push(bwtintv_t, *curr, ik);
            break; // always terminate extension at an ambiguous base; in this case, i<len always stands
        }
    }

    if (i == len) { // since we do the backwards search as the forwards search
        kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
    }
    bwt_reverse_intvs(curr); // s.t. smaller intervals (i.e. longer matches) visited first
    ret = curr->a[0].info; // this will be the returned value
    swap = curr; curr = prev; prev = swap;

    for (i = x - 1; i >= -1; --i) { // backward search for MEMs
        c = i < 0? -1 : q[i] < 4? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base

        for (j = 0, curr->n = 0; (unsigned) j < prev->n ; ++j) {
            bwtintv_t *p = &prev->a[j];

            if (c >= 0 && ik.x[2] >= max_intv)
                bwt_extend(bwt, p, ok, 1, c); // is_back == true

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
// this copy of bwt_extend works with bwt_smem1a
// Called from memchain.c
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
            bwt_extend(bwtc, &ik, ok, 0, c);
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

void bwt_dump_bwt(const char *fn, const bwt_t *bwt, ubyte_t I) {
    FILE *fp;
    fp = xopen(fn, "wb");
    err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
    err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    //second argument is the number of bytes; third argument is the number of elements
    err_fwrite(bwt->bwt, 4, bwt->bwt_size, fp);
    err_fflush(fp);
    err_fclose(fp);
}

void bwt_dump_bwt_new(const char *fn, const bwt_t *bwt, ubyte_t I) {
    FILE *fp;
    fp = xopen(fn, "wb");
    err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
    err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    //second argument is the number of bytes; third argument is the number of elements
    err_fwrite(bwt->bwt_new, 8, bwt->bwt_size, fp);
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

// bwt_restore_bwt and bwt_restore_bwt2 were changed so that they restore the correct data structure (bwt->bwt).
bwt_t *bwt_restore_bwt(const char *fn) {
    bwt_t *bwt;
    FILE *fp;

    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    fp = xopen(fn, "rb");
    err_fseek(fp, 0, SEEK_END);
    bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
    bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 8);
    err_fseek(fp, 0, SEEK_SET);
    err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
    err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    fread_fix(fp, bwt->bwt_size<<2, bwt->bwt); //change bwt->bwt_size<<2 with seq_len/8
    bwt->seq_len = bwt->L2[4];
    err_fclose(fp);

    return bwt;
}

bwt_t *bwt_restore_bwt_new(const char *fn) {
    bwt_t *bwt;
    FILE *fp;

    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    fp = xopen(fn, "rb");
    err_fseek(fp, 0, SEEK_END);
    bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
    bwt->bwt_new = (uint64_t*)calloc(bwt->bwt_size, 8);
    err_fseek(fp, 0, SEEK_SET);
    err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
    err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    fread_fix(fp, bwt->bwt_size<<2, bwt->bwt_new); //change bwt->bwt_size<<2 with seq_len/8
    bwt->seq_len = bwt->L2[4];
    err_fclose(fp);

    return bwt;
}
// whats the difference between bwt_restore_bwt and bwt_restore_bwt2?
// bwt_restore_bwt2 is used because bwt is already connected to idx
void bwt_restore_bwt2(const char *fn, bwt_t *bwt) {
    FILE *fp;
    memset(bwt, 0, sizeof(bwt_t));
    /* bwt = (bwt_t*)calloc(1, sizeof(bwt_t)); */
    fp = xopen(fn, "rb");
    err_fseek(fp, 0, SEEK_END);
    bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
    bwt->bwt_new = (uint64_t*)calloc(bwt->bwt_size, 8);
    err_fseek(fp, 0, SEEK_SET);
    err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
    err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    fread_fix(fp, bwt->bwt_size<<2, bwt->bwt_new);
    bwt->seq_len = bwt->L2[4];

    err_fclose(fp);
}

void bwt_destroy(bwt_t *bwt) {
    if (bwt == 0) return;
    if (bwt->sa) free(bwt->sa);
    if (bwt->bwt_new) free(bwt->bwt_new);
    if (bwt->bwt) free(bwt->bwt);
    free(bwt);
}
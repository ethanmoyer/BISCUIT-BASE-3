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
bwtint_t builtin_popcountll(uint64_t seq0, uint64_t seq1, int c, bwtint_t k, uint8_t parent) {
    if (k == 0) // How can we remove this...
        return 0;
    if (k > 64)
        k = 64;
    int shift = 64 - k;
    if (parent) {
        switch(c) { // G>A
            case 2:
                return 0;
            case 0:
                return __builtin_popcountll(~seq1 >> shift); //A
            case 1:
                return __builtin_popcountll(seq0 >> shift & seq1 >> shift); //C
        }
        return __builtin_popcountll(~seq0 >> shift); //T
    } else {
        switch(c) { // C>T
            case 1:
                return 0;
            case 0:
                return __builtin_popcountll(seq0 >> shift); //A
            case 2:
                return __builtin_popcountll(~seq0 >> shift & ~seq1 >> shift); //G
        }
        return __builtin_popcountll(seq1 >> shift); //T
    }
    return 0;
}

bwtint_t bwt_occ_new_index(const bwt_t *bwt, bwtint_t k, int c, uint8_t parent) {
    // bwt starts indexing at 0, so calculations involving k are handled accordingly.
    if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
    if (k == (bwtint_t)(-1)) return 0;

    k -= (k >= bwt->primary);
    uint32_t index = k/128 * 8;
    k++;

    bwtint_t count = 0;
    if (k / 128)
        count = bwt->bwt_new[(k / 128 - 1) * 8 + c];
    if (k & 127 == 0)
        return count;

    count += builtin_popcountll(bwt->bwt_new[index + 4], bwt->bwt_new[index + 6], c, k & 127, parent);
    if ((k & 127) > 64)
        count += builtin_popcountll(bwt->bwt_new[index + 5], bwt->bwt_new[index + 7], c, k & 63, parent);

    return count;
}

// compute inverse CSA
static inline bwtint_t bwt_invPsi(bwt_t *bwt, bwtint_t k, uint8_t parent) {
    bwtint_t x = k - (k > bwt->primary);
    int t = x;
    //fprintf(stderr, "k: %d\n", k);
    //x = nucAtBWTinv_new(bwt, x, !parent);
    //fprintf(stderr, "x: %d\n", x);
    int i = 8 * (t/128) + 4 + ((t & 127) >= 64);
    //fprintf(stderr, "index1: %d index2: %d\n", i, i + 2);

    //bwtint_t bwt0 = bwt->bwt_new[i];
    //bwtint_t bwt1 = bwt->bwt_new[i + 2];

    //fprintf(stderr, "bwt0 \t %llu \t %s\n", bwt0, "1111111111111100000011001100111111111100000001101110000000011000");
    //fprintf(stderr, "bwt1 \t %llu \t %s\n", bwt1, "0001100111111000000011110000100000000011110000000000111100100001");

    // 111001100000000000000000 0011011011111100001111110011000000000000
    // 000110011111100000001111 0000100000000011110000000000111100100001
    int w = (~(t)&0x3f);
    //fprintf(stderr, "w: %d\n", w);
    int y = (bwt->bwt_new[i]>>w)&1 ^ 1;
    int z = (bwt->bwt_new[i + 2]>>w)&1 & 1;

    //fprintf(stderr, "y: %d\n", y);
    //fprintf(stderr, "z: %d\n", z);
    int q = (y << 1) | z;
    //fprintf(stderr, "q: %d\n\n", q);
    //if (x != q) {
    //    fprintf(stderr, "ERROR\n");
    //    exit(0);
    //}

    x = q;
    x = bwt->L2[x] + bwt_occ_new_index(bwt, k, x, !parent);
    return k == bwt->primary ? 0 : x;
}


void bwt_cal_sa(bwt_t *bwt, int intv, uint8_t parent)
{
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
            bwt->sa[isa / intv] = sa;
        }
        --sa;
        isa = bwt_invPsi(bwt, isa, parent);
    }

    if (isa % intv == 0) bwt->sa[isa/intv] = sa;
    bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

bwtint_t bwt_sa(bwt_t *bwt, bwtint_t k, uint8_t parent) {
    bwtint_t sa = 0, mask = bwt->sa_intv - 1;

    while (k & mask) { // != 0
        ++sa;
        k = bwt_invPsi(bwt, k, parent);
    }

    /* without setting bwt->sa[0] = -1, the following line should be
       changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
    return sa + bwt->sa[k/bwt->sa_intv];
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
void bwt_extend_debug(const bwt_t *bwt, bwtintv_t *ik, bwtintv_t *ok, int is_back, int c, uint8_t parent) {
    bwtint_t tk[4] = {0};
    bwtint_t tl[4] = {0};
    uint8_t dir = is_back ^ parent;
    uint8_t _is_back = !is_back;
    for (int i = 3; i != c - 1; i--) {
        tk[i] = bwt_occ_new_index(bwt, ik->x[_is_back] - 1, i, dir);
        tl[i] = bwt_occ_new_index(bwt, ik->x[_is_back] + ik->x[2] - 1, i, dir);
        ok[i].x[_is_back] = bwt->L2[i] + 1 + tk[i];
        ok[i].x[2] = tl[i] - tk[i];
        if (i == 3)
            ok[3].x[is_back] = ik->x[is_back] + (
                    ik->x[_is_back] <= bwt->primary && ik->x[_is_back] + ik->x[2] - 1 >= bwt->primary);
        if (i != 0)
            ok[i - 1].x[is_back] = ok[i].x[is_back] + ok[i].x[2];
    }
}


int bwt_smem1a(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2], uint8_t parent) {

    printf(stderr, "[bwt_snem1a] Start\n");
    double t_real;
    t_real = realtime();

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
        if (ik.x[2] < max_intv) { // an interval small enough, currently won't come here
            kv_push(bwtintv_t, *curr, ik);
            break;
        } else if (q[i] < 4) { // an A/C/G/T base
            c = 3 - q[i]; // complement of q[i]
            bwt_extend_debug(bwtc, &ik, ok, 0, c, parent);
            if (ok[c].x[2] != ik.x[2]) { // change of the interval size
                kv_push(bwtintv_t, *curr, ik);

                // no more matches
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
            if (c >= 0 && ik.x[2] >= max_intv) bwt_extend_debug(bwt, p, ok, 1, c, parent);
            //fprintf(stderr, "ok[c].x[2]: %llu\n", ok[c].x[2]);
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

    fprintf(stderr, "\n[bwt_snem1a] Real time: %.3f sec\n", realtime()-t_real);

    return ret;
}
// this copy of bwt_extend works with bwt_smem1a
// Called from memchain.c
int bwt_smem1(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2], uint8_t parent) {
    return bwt_smem1a(bwt, bwtc, len, q, x, min_intv, 0, mem, tmpvec, parent);
}

int bwt_seed_strategy1(const bwt_t *bwt, const bwt_t *bwtc, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem, uint8_t parent) {
    int i, c;
    bwtintv_t ik, ok[4];

    memset(mem, 0, sizeof(bwtintv_t));
    if (q[x] > 3) return x + 1;
    bwt_set_intv(bwt, bwtc, q[x], ik); // the initial interval of a single base
    for (i = x + 1; i < len; ++i) { // forward search
        if (q[i] < 4) { // an A/C/G/T base
            c = 3 - q[i]; // complement of q[i]
            bwt_extend_debug(bwtc, &ik, ok, 0, c, parent);
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
    //second argument is the number of bytes; third argument is the number of elements
    err_fwrite(bwt->bwt, 4, bwt->bwt_size, fp);
    err_fflush(fp);
    err_fclose(fp);
}

void bwt_dump_bwt_new(const char *fn, const bwt_t *bwt) {
    FILE *fp;
    fp = xopen(fn, "wb");
    err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
    err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    err_fwrite(bwt->bwt_new, 8, (bwt->seq_len/128 + 1) * 16, fp);
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
    bwt->seq_len = err_ftell(fp) * 2;
    bwt->bwt_new = (uint64_t*)calloc((bwt->seq_len/128 + 1) * 16, 8);
    err_fseek(fp, 0, SEEK_SET);
    err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
    err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    bwt->seq_len = bwt->L2[4];
    fread_fix(fp, bwt->bwt_size<<2, bwt->bwt_new); //change bwt->bwt_size<<2 with seq_len/8

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
    bwt->seq_len = err_ftell(fp) * 2;
    bwt->bwt_new = (uint64_t*)calloc((bwt->seq_len/128 + 1) * 16, 8);
    err_fseek(fp, 0, SEEK_SET);
    err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
    err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
    fread_fix(fp, bwt->bwt_size<<2, bwt->bwt_new);
    bwt->seq_len = bwt->L2[4];
    //fprintf(stderr, "bwt->bwt_size: %llu bwt->seq_len: %llu\n", bwt->bwt_size, bwt->seq_len);


    err_fclose(fp);
}

void bwt_destroy(bwt_t *bwt) {
    if (bwt == 0) return;
    if (bwt->sa) free(bwt->sa);
    if (bwt->bwt_new) free(bwt->bwt_new);
    if (bwt->bwt) free(bwt->bwt);
    free(bwt);
}
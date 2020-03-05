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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include <limits.h>
#include <math.h>
#include "bntseq.h"
#include "bwt.h"
#include "utils.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int is_bwt(ubyte_t *T, int n);

int64_t bwa_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is)
{
    bwt_t *bwt;
    ubyte_t *buf, *buf2;
    uint32_t i, pac_size;         /* WZ from int */
    FILE *fp;

    // initialization
    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    bwt->seq_len = bwa_seq_len(fn_pac);
    bwt->bwt_size = (bwt->seq_len + 15) >> 4;
    fp = xopen(fn_pac, "rb");

    // prepare sequence
    pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
    buf2 = (ubyte_t*)calloc(pac_size, 1);
    err_fread_noeof(buf2, 1, pac_size, fp);
    err_fclose(fp);
    memset(bwt->L2, 0, 5 * 4);
    buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
    for (i = 0; i < bwt->seq_len; ++i) {
        buf[i] = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
        ++bwt->L2[1+buf[i]];
    }
    for (i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];
    free(buf2);

    // Burrows-Wheeler Transform
    if (use_is) {
        bwt->primary = is_bwt(buf, bwt->seq_len);
    } else {
#ifdef _DIVBWT
        bwt->primary = divbwt(buf, buf, 0, bwt->seq_len);
#else
        err_fatal_simple("libdivsufsort is not compiled in.");
#endif
    }

    bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
    for (i = 0; i < bwt->seq_len; ++i)
        bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);

    free(buf);

    return bwt;
}

int bwa_pac2bwt(int argc, char *argv[]) // the "pac2bwt" command; IMPORTANT: bwt generated at this step CANNOT be used with BWA. bwtupdate is required!
{
	bwt_t *bwt;
	int c, use_is = 1;
	while ((c = getopt(argc, argv, "d")) >= 0) {
		switch (c) {
		case 'd': use_is = 0; break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa pac2bwt [-d] <in.pac> <out.bwt>\n");
		return 1;
	}
	bwt = bwt_pac2bwt(argv[optind], use_is);
	bwt_dump_bwt(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

//Where the occurrence is generated every OCC_INTERVAL
// This function generates the counts after each string of 128 bases. Using the __builtin_popcountll function, strings
// of 64 nucleotides are counting with one operation. Depending on the index int that is passed, either C's (index = 0)
// or G's (index = 1) will be ignored. The accessing of the bwt in bwt->bwt_new falls in line with the eight block
// data structure.
// The output after each 128 block is simply used for developmental purposes and debugging of the counts.
void bwt_bwtupdate_core(bwt_t *bwt, int index) {
    bwtint_t n_occ;
    n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
    bwt->bwt_new = (uint64_t*) calloc((bwt->seq_len / 128 + 1) * 16 , 8);
    bwt->bwt_size += n_occ * sizeof(bwtint_t);
    uint64_t j = 1UL;
    uint64_t x = 0;
    ubyte_t buf = 0;

    for (bwtint_t i = 0; i < bwt->seq_len + 16; ++i) {
        buf = 0;
        buf = (bwt->bwt[i/16] >> (30 - 2 * i % 32)) & 3;
        x = j << (63 - (i % 64));
        if (buf <= 1) {
            bwt->bwt_new[i/128 * 8 + (i % 128 < 64 ? 4 : 5)] ^= x;
        }
        if (!(buf % 2 == 0)) {
            bwt->bwt_new[i/128 * 8 + (i % 128 < 64 ? 6 : 7)] ^= x;
        }
    }

    bwtint_t n = bwt->seq_len / 128 + 1; //1400 --> 700
    bwtint_t k = 8;
    for (bwtint_t i = 0; i < n - 1; i++) {

        //A
        bwt->bwt_new[i * k + 0] =
            __builtin_popcountll((bwt->bwt_new[i * k + 4] ^ bwt->bwt_new[i * k + 6]) & bwt->bwt_new[i * k + 4]) +
            __builtin_popcountll((bwt->bwt_new[i * k + 5] ^ bwt->bwt_new[i * k + 7]) & bwt->bwt_new[i * k + 5]);

        //T
        bwt->bwt_new[i * k + 3] =
            __builtin_popcountll(bwt->bwt_new[i * k + 6] & ~bwt->bwt_new[i * k + 4]) +
            __builtin_popcountll(bwt->bwt_new[i * k + 7] & ~bwt->bwt_new[i * k + 5]);

        if (index == 0) {
            //G
            bwt->bwt_new[i * k + 2] =
                    __builtin_popcountll(~(bwt->bwt_new[i * k + 4] ^
                                           bwt->bwt_new[i * k + 6])) +
                    __builtin_popcountll(~(bwt->bwt_new[i * k + 5] ^
                                           bwt->bwt_new[i * k + 7]));
        } else {
            //C
            bwt->bwt_new[i * k + 1] =
                    __builtin_popcountll(bwt->bwt_new[i * k + 4] &
                                         bwt->bwt_new[i * k + 6]) +
                    __builtin_popcountll(bwt->bwt_new[i * k + 5] &
                                         bwt->bwt_new[i * k + 7]);
        }

        if (i != 0) {
            bwt->bwt_new[i * k + 0] += bwt->bwt_new[(i - 1) * k + 0];
            bwt->bwt_new[i * k + 2] += bwt->bwt_new[(i - 1) * k + 2];
            bwt->bwt_new[i * k + 3] += bwt->bwt_new[(i - 1) * k + 3];
            bwt->bwt_new[i * k + 1] += bwt->bwt_new[(i - 1) * k + 1];
        }

    }
}

// There's going to be a problem here if bwtupdate is run because bwt_bwtupdate_core needs to be run twice in order to
// have both G>A and C>T
// Asking wanding about how he wants to handle this functions and functions alike. We only store enough data for one of
// the indexes and not two. How was it handled before.
int bwa_bwtupdate(int argc, char *argv[]) // the "bwtupdate" command
{
    bwt_t *bwt;
    if (argc < 2) {
        fprintf(stderr, "Usage: bwa bwtupdate <the.bwt>\n");
        return 1;
    }
    bwt = bwt_restore_bwt(argv[1]);
    bwt_bwtupdate_core(bwt, 0);
    bwt_dump_bwt(argv[1], bwt);
    bwt_destroy(bwt);
    return 0;
}


int bwa_bwt2sa(int argc, char *argv[]) // the "bwt2sa" command
{
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
		case 'i': sa_intv = atoi(optarg); break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa bwt2sa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_sa(bwt, sa_intv, 1);
	bwt_dump_sa(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

int main_biscuit_index(int argc, char *argv[]) {

    extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

    char *prefix = 0, *str, *str2, *str3;
    int c, algo_type = 0, is_64 = 0;
    clock_t t;
    int64_t l_pac;

    while ((c = getopt(argc, argv, "6a:p:")) >= 0) {
        switch (c) {
        case 'a': // if -a is not set, algo_type will be determined later
            if (strcmp(optarg, "div") == 0) algo_type = 1;
            else if (strcmp(optarg, "bwtsw") == 0) algo_type = 2;
            else if (strcmp(optarg, "is") == 0) algo_type = 3;
            else err_fatal(__func__, "unknown algorithm: '%s'.", optarg);
            break;
        case 'p': prefix = strdup(optarg); break;
        case '6': is_64 = 1; break;
        default: return 1;
        }
    }

    if (optind + 1 > argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   biscuit index [-a bwtsw|is] [-c] <in.fasta>\n\n");
        fprintf(stderr, "Options: -a STR    BWT construction algorithm: bwtsw or is [auto]\n");
        fprintf(stderr, "         -p STR    prefix of the index [same as fasta name]\n");
        fprintf(stderr, "         -6        index files named as <in.fasta>.64.* instead of <in.fasta>.* \n");
        fprintf(stderr, "\n");
        fprintf(stderr,	"Warning: `-a bwtsw' does not work for short genomes, while `-a is' and\n");
        fprintf(stderr, "         `-a div' do not work not for long genomes. Please choose `-a'\n");
        fprintf(stderr, "         according to the length of the genome.\n\n");
        return 1;
    }
    if (prefix == 0) {
        prefix = malloc(strlen(argv[optind]) + 4);
        strcpy(prefix, argv[optind]);
        if (is_64) strcat(prefix, ".64");
    }
    str  = (char*)calloc(strlen(prefix) + 50, 1);
    str2 = (char*)calloc(strlen(prefix) + 50, 1);
    str3 = (char*)calloc(strlen(prefix) + 50, 1);

    {   /* nucleotide indexing */
        /* generates .ct.pac, .ct.ann and .ct.amb */
        gzFile fp = xzopen(argv[optind], "r");
        t = clock();
        fprintf(stderr, "[%s] Pack bisulfite FASTA... ", __func__);
        l_pac = bis_bns_fasta2bntseq(fp, prefix, 1); /* parent strand */
        l_pac = bis_bns_fasta2bntseq(fp, prefix, 0); /* daughter strand */

        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
        err_gzclose(fp);
    }
    // set the algorithm for generating BWT
    if (algo_type == 0) algo_type = l_pac > 50000000? 2 : 3; // set the algorithm for generating BWT
    {
        t = clock();
        fprintf(stderr, "[%s] Construct BWT for the parent strands...\n", __func__);
        strcpy(str, prefix); strcat(str, ".par.pac");
        strcpy(str2, prefix); strcat(str2, ".par.bwt");
        if (algo_type == 2) bwt_bwtgen(str, str2);
        else if (algo_type == 1 || algo_type == 3) {
            bwt_t *bwt;
            bwt = bwt_pac2bwt(str, algo_type == 3);
            bwt_dump_bwt(str2, bwt);
            bwt_destroy(bwt);

        }
        fprintf(stderr, "[%s] %.2f seconds elapse.\n", __func__, (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    {
        t = clock();
        fprintf(stderr, "[%s] Construct BWT for the daughter strands...\n", __func__);
        strcpy(str, prefix); strcat(str, ".dau.pac");
        strcpy(str2, prefix); strcat(str2, ".dau.bwt");
        if (algo_type == 2) bwt_bwtgen(str, str2);
        else if (algo_type == 1 || algo_type == 3) {
            bwt_t *bwt;
            bwt = bwt_pac2bwt(str, algo_type == 3);
            bwt_dump_bwt(str2, bwt);

            bwt_destroy(bwt);
        }
        fprintf(stderr, "[%s] %.2f seconds elapse.\n", __func__, (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    {
        bwt_t *bwt;
        strcpy(str, prefix); strcat(str, ".par.bwt");
        t = clock();
        fprintf(stderr, "[%s] Update parent BWT... \n", __func__);
        bwt = bwt_restore_bwt(str);
        bwt_bwtupdate_core(bwt, 0);
        bwt_dump_bwt_new(str, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "[%s] %.2f sec\n", __func__, (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    {
        bwt_t *bwt;
        strcpy(str, prefix); strcat(str, ".dau.bwt");
        t = clock();
        fprintf(stderr, "[%s] Update daughter BWT... \n", __func__);
        bwt = bwt_restore_bwt(str);
        bwt_bwtupdate_core(bwt, 1);
        bwt_dump_bwt_new(str, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "[%s] %.2f sec\n", __func__, (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    {
        gzFile fp = xzopen(argv[optind], "r");
        t = clock();
        fprintf(stderr, "[%s] Pack forward-only FASTA... ", __func__);
        l_pac = dump_forward_pac(fp, prefix);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
        err_gzclose(fp);
        strcpy(str, prefix); strcat(str, ".par.pac");
        unlink(str);
        strcpy(str, prefix); strcat(str, ".dau.pac");
        unlink(str);
    }
    {
        bwt_t *bwt;
        strcpy(str, prefix); strcat(str, ".par.bwt");
        strcpy(str3, prefix); strcat(str3, ".par.sa");
        t = clock();
        fprintf(stderr, "[%s] Construct parent SA from BWT and Occ... ", __func__);
        bwt = bwt_restore_bwt_new(str);
        bwt_cal_sa(bwt, 32, 1);
        bwt_dump_sa(str3, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    {
        bwt_t *bwt;
        strcpy(str, prefix); strcat(str, ".dau.bwt");
        strcpy(str3, prefix); strcat(str3, ".dau.sa");
        t = clock();
        fprintf(stderr, "[%s] Construct daughter SA from BWT and Occ... ", __func__);
        bwt = bwt_restore_bwt_new(str);
        bwt_cal_sa(bwt, 32, 0);
        bwt_dump_sa(str3, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    free(str3); free(str2); free(str); free(prefix);
    return 0;
}

/* lower case repetitive region */

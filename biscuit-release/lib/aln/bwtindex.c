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

unsigned long long x;
unsigned long long j = 1UL;

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

//return array of two bwt_t
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
	fprintf(stderr, "BWT size: %llu\n", bwt->bwt_size);
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
	buf2 = (ubyte_t*)calloc(pac_size, 1);
	err_fread_noeof(buf2, 1, pac_size, fp);
	err_fclose(fp);
	memset(bwt->L2, 0, 5 * 4);
	buf = (ubyte_t*)calloc((bwt->seq_len + 1)/2, 1);
	for (i = 0; i < bwt->seq_len/2; ++i) {
		buf[i] = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
		++bwt->L2[1+(buf2[i>>2] >> ((3 - (i&3)) << 1) & 3)];
	}
	for (i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];

	free(buf2);
    //for (i = 0; i < bwt->seq_len/2; i++)
    //    fprintf(stderr, "buf: %d\n", buf[i]);
	// Burrows-Wheeler Transform
	if (use_is) {
		bwt->primary = is_bwt(buf, bwt->seq_len/2);
    fprintf(stderr, "\n");
	//for (i = 0; i < bwt->seq_len/2; i++)
	//	fprintf(stderr, "buf: %d\n", buf[i]);

	} else {
#ifdef _DIVBWT
		bwt->primary = divbwt(buf, buf, 0, bwt->seq_len);

#else
		err_fatal_simple("libdivsufsort is not compiled in.");
#endif
	}

	//allocate size of vector
    bwt->bwt_all.size = ceil((double) bwt->seq_len/128);

    // intialize space of vector
    bwt->bwt_all.bwt0 = (uint64_t*)calloc(bwt->bwt_all.size, 8);
    bwt->bwt_all.bwt1 = (uint64_t*)calloc(bwt->bwt_all.size, 8);


    for (i = 0; i < bwt->bwt_all.size; i++) {
        bwt->bwt_all.bwt0[i] = ULLONG_MAX;
        bwt->bwt_all.bwt1[i] = ULLONG_MAX;
    }

    for (i = 0; i < bwt->seq_len/2; ++i) {
        x = j << (63 - (i % 64));
        if (buf[i] > 1) {
            bwt->bwt_all.bwt0[(int)ceil(i/64)] ^= x;
        }

        if (buf[i] % 2 == 0) {
            bwt->bwt_all.bwt1[(int)ceil(i/64)] ^= x;
        }
    }

    bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
    for (i = 0; i < bwt->seq_len; ++i)
    	bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);

    free(buf);
    return bwt;

}

uint64_t* newIndex(uint64_t s, int length) {
    uint64_t *s_ = (uint64_t*)calloc(2, 8);
    s_[0] = ULLONG_MAX;
    s_[1] = ULLONG_MAX;
    for (int i = 0; i < length; i++) {
        x = j << (length - 1 - i);
        //fprintf(stderr, "x: %llu\n", x);
        //fprintf(stderr, "i/2: %d s >> i & 1: %d\n", i/2, (s & x) > 1);
        //fprintf(stderr,"x: %llu ((s & x) > 1: %d\n", x, ((s & x) >= 1));

        if (i % 2 == 0) {
            if (((s & x) >= 1)) {
                s_[0] ^= j << (63 - i/2); //this needs to change
                //fprintf(stderr, "i: %d\n", i);
            }
        } else {
            if (!((s & x) >= 1)) {
                s_[1] ^= j << (63 - i/2);
            }
        }
    }
    return s_;
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
//(I think that’s 128 bases) and write to the new bwt (which is buf)
void bwt_bwtupdate_core(bwt_t *bwt)
{
    bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	//fprintf(stderr, "n_occ: %llu\n", n_occ);
	//exit(0);

	bwt->bwt_size += n_occ * sizeof(bwtint_t); // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i) {
		if (i % OCC_INTERVAL == 0) {
			memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += sizeof(bwtint_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0) buf[k++] = bwt->bwt[i/16]; // 16 == sizeof(uint32_t)/2
		++c[bwt_B00(bwt, i)];
	}
	// the last element
	memcpy(buf + k, c, sizeof(bwtint_t) * 4);
	xassert(k + sizeof(bwtint_t) == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt); bwt->bwt = buf;

    bwt->bwt_occ_matrix0.rows = floor(bwt->seq_len / 256);
    fprintf(stderr, "\nCheckpoints!\n");

    bwt->bwt_occ_matrix0.cols = 3 ;

    bwt->bwt_occ_matrix0.occurrences = (uint64_t*)
    calloc(bwt->bwt_occ_matrix0.rows * bwt->bwt_occ_matrix0.cols, 8);

    for (i = 0; i < bwt->bwt_all.size; i++)  {
        fprintf(stderr, "bwt0: %llu \n", bwt->bwt_all.bwt0[i]);
        fprintf(stderr, "bwt1: %llu \n", bwt->bwt_all.bwt1[i]);
    }

    i = 0;
    j = 1;
	for (k = 0; k < bwt->bwt_occ_matrix0.rows; ) {

        //G
        bwt->bwt_occ_matrix0.occurrences[k * 3 + 0] =
            __builtin_popcountll(bwt->bwt_all.bwt0[i]) +
            __builtin_popcountll(bwt->bwt_all.bwt0[j]);
        //T
        bwt->bwt_occ_matrix0.occurrences[k * 3 + 1] =
            __builtin_popcountll(bwt->bwt_all.bwt1[i]) +
            __builtin_popcountll(bwt->bwt_all.bwt1[j]);
        //A
        bwt->bwt_occ_matrix0.occurrences[k * 3 + 2] =
            __builtin_popcountll(~(bwt->bwt_all.bwt1[i] ^
                    bwt->bwt_all.bwt0[i])) +
            __builtin_popcountll(~(bwt->bwt_all.bwt1[j] ^
                    bwt->bwt_all.bwt0[j]));

        if (k != 0) {
            bwt->bwt_occ_matrix0.occurrences[k * 3 + 0] +=
                bwt->bwt_occ_matrix0.occurrences[(k - 1) * 3 + 0];
            bwt->bwt_occ_matrix0.occurrences[k * 3 + 1] +=
                 bwt->bwt_occ_matrix0.occurrences[(k - 1) * 3 + 1];
            bwt->bwt_occ_matrix0.occurrences[k * 3 + 2] +=
                 bwt->bwt_occ_matrix0.occurrences[(k - 1) * 3 + 2];
        }

        fprintf(stderr, "k: %llu occurances: A: %llu T: %llu G: %llu \n", k,
        bwt->bwt_occ_matrix0.occurrences[k * 3 + 0],
        bwt->bwt_occ_matrix0.occurrences[k * 3 + 1],
        bwt->bwt_occ_matrix0.occurrences[k * 3 + 2]);

        k++;
        i+=2;
        j+=2;
	}

    //bwt->occurance = o;
    fprintf(stderr, "\n");

    free(bwt->bwt_occ_matrix0.occurrences);

}

int bwa_bwtupdate(int argc, char *argv[]) // the "bwtupdate" command
{
	bwt_t *bwt;
	if (argc < 2) {
		fprintf(stderr, "Usage: bwa bwtupdate <the.bwt>\n");
		return 1;
	}
	bwt = bwt_restore_bwt(argv[1]);
	bwt_bwtupdate_core(bwt);
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
	bwt_cal_sa(bwt, sa_intv);
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

  { /* nucleotide indexing */
    /* generates .ct.pac, .ct.ann and .ct.amb */
    gzFile fp = xzopen(argv[optind], "r");
    t = clock();
    fprintf(stderr, "[%s] Pack bisulfite FASTA... ", __func__);
    l_pac = bis_bns_fasta2bntseq(fp, prefix, 1); /* parent strand */
    l_pac = bis_bns_fasta2bntseq(fp, prefix, 0); /* daughter strand */

    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    err_gzclose(fp);
  }
  if (algo_type == 0) algo_type = l_pac > 50000000? 2 : 3; // set the algorithm for generating BWT
  {
    t = clock();
    fprintf(stderr, "[%s] Construct BWT for the parent strands...\n", __func__);
    strcpy(str, prefix); strcat(str, ".par.pac");
    strcpy(str2, prefix); strcat(str2, ".par.bwt");
    if (algo_type == 2) bwt_bwtgen(str, str2);
    else if (algo_type == 1 || algo_type == 3) {
        bwt_t *bwt;

        // this is the subsequence that we're searching with old index
        // 00 = A, 10 = G, 11 = T
        uint64_t s = 2814;
        // the length var below is to counteract the missing data from the 00 of As
        int s_len = 16;

        //store subseq with new index
        uint64_t *subseq = newIndex(s, s_len);

        fprintf(stderr, "subseq0: %llu subseq1: %llu\n", subseq[0], subseq[1]);

        //generate bwt
        bwt = bwt_pac2bwt(str, algo_type == 3);
        //generates old occurrences array
        bwt_gen_cnt_table(bwt);
        //generates new occurrences array
        bwt_bwtupdate_core(bwt);

        bwtint_t k = 165;

        //A = 0
        //T = 1
        //G = 2

        bwtint_t cntk[4];

        //this changes bwt.. why?
        bwt_occ4_new_index(bwt, bwt->seq_len/2, cntk);

        //suffix array
        bwt_cal_sa(bwt, 128);

        fprintf(stderr, "primary %llu\n", bwt->primary);

        k = 128;
        bwtint_t t = bwt_sa(bwt, 128);

        fprintf(stderr, "t %llu\n", t);

        // bwtc is the bwt on the comp. strand --> bi directional search
        // q -> base sequence; base-3?
        // ik -> input interval for iteration
        // ok -> output interval for iteration

        strcpy(str, prefix); strcat(str, ".dau.pac");
/*
        bwt_t *bwtc;
        bwtc = bwt_pac2bwt(str, algo_type == 3);
*/
        //A..., C..., G..., T.... and $
        uint8_t q = 1; //AAGG? nope 0 through 4
        bwtintv_t ik, ok[3];

        fprintf(stderr, "L2[0]: %llu\n", bwt->L2[0]);
        fprintf(stderr, "L2[1]: %llu\n", bwt->L2[1]);
        fprintf(stderr, "L2[3]: %llu\n", bwt->L2[3]);
        fprintf(stderr, "L2[4]: %llu\n", bwt->L2[4]);

        //L2[1] - L2[0] = 33      //#A
        //L2[2] - L2[1] = 0       //#C
        //L2[3] - L2[2] = 96      //#G
        //L2[3] - L2[4] = 128     //#T

        /**
         * x[0] - forward index location;
         * x[1] - reverse complement index location;
         * x[2] - occurrence number (number of substrings);
         * info>>32 - beg;
         * (uint32_t) info - end;
         */

         /*     0 = A
          *     2 = G
          *     3 = T
          */
        //bwt_set_intv(bwt, bwt, 2, ik);
        //bwt_set_intv(bwt, bwt, 3, ik); // the initial interval of a single base
        fprintf(stderr, "index: %llu occ: %llu\n", ik.x[0], ik.x[2]);

        //const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back
        bwt_extend(bwt, &ik, ok, 1, subseq, s_len/2);
exit(0);
        //bwt_extend(bwt, &ik, ok, 2);


        //bwt_extend(bwt, &ik, ok, 0);

        bwt_dump_bwt(str2, bwt);
        bwt_destroy(bwt);
        exit(0);
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
    bwt_bwtupdate_core(bwt);
    bwt_dump_bwt(str, bwt);
    bwt_destroy(bwt);
    fprintf(stderr, "[%s] %.2f sec\n", __func__, (float)(clock() - t) / CLOCKS_PER_SEC);
  }
  {
    bwt_t *bwt;
    strcpy(str, prefix); strcat(str, ".dau.bwt");
    t = clock();
    fprintf(stderr, "[%s] Update daughter BWT... \n", __func__);
    bwt = bwt_restore_bwt(str);
    bwt_bwtupdate_core(bwt);
    bwt_dump_bwt(str, bwt);
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
    bwt = bwt_restore_bwt(str);
    bwt_cal_sa(bwt, 32);
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
    bwt = bwt_restore_bwt(str);
    bwt_cal_sa(bwt, 32);
    bwt_dump_sa(str3, bwt);
    bwt_destroy(bwt);
    fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  }
  free(str3); free(str2); free(str); free(prefix);
  return 0;
}

/* lower case repetitive region */

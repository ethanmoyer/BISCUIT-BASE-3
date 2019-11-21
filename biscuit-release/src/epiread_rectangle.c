/* Rectangularize an epiread output */

#include "wztsv.h"
#include "wvec.h"
#include "kstring.h"
#include "refcache.h"

int refcache_next_cg(refcache_t *rc, int pos) {
   for (;; pos++) {
      if (toupper(refcache_getbase_auto(rc, NULL, pos)) == 'C' &&
          toupper(refcache_getbase_auto(rc, NULL, pos+1)) == 'G')
         return pos;
   }
}

static int usage() {
   fprintf(stderr, "\n");
   fprintf(stderr, "Convert epiread into a rectanglular matrix.\n");
   fprintf(stderr, "Usage: biscuit rectangle [options] [ref.fa] [.epiread]\n");
   fprintf(stderr, "Options:\n\n");
   fprintf(stderr, "     -o        output file [stdout]\n");
   fprintf(stderr, "     -h        this help.\n");
   fprintf(stderr, "\n");
   return 1;
}

typedef struct read_t {
   kstring_t seq;
   kstring_t other; // other fields of the line
} read_t;

DEFINE_VECTOR(read_v, read_t)

int main_rectangle(int argc, char *argv[]) {

   char *out_fn = 0;
   
   if (argc<2) return usage();
   int c;
   while ((c = getopt(argc, argv, "o:"))>=0) {
      switch (c) {
      case 'o': out_fn = optarg; break;
      case 'h': return usage();
      default:
         fprintf(stderr, "Unrecognized command\n");
      }
   }

   char *ref_fn = argv[optind++];
   char *epiread_fn = argv[optind++];
   refcache_t *rc = init_refcache(ref_fn, 1000, 1000);

   tsv_t *tsv = tsv_open(epiread_fn);
   int region_beg = 0, region_width = -1;
   read_v *reads = init_read_v(1000);
   char *chrm = NULL;
   while (tsv_read(tsv)) {
      if (tsv_is_blankline(tsv)) continue;
      char *f = tsv_field(tsv, 4);
      if (f[0] == '.') {
         read_t *r = next_ref_read_v(reads);
         r->seq = (const kstring_t) {0};
         r->other = (const kstring_t) {0};
         kputs(tsv->line, &r->other);
         continue;
      }
      
      int read_beg = atoi(f);
      if (!region_beg) region_beg = read_beg;

      if (chrm == NULL) {
         chrm = calloc(strlen(tsv_field(tsv, 0))+1, sizeof(char));
         strcpy(chrm, tsv_field(tsv, 0));
         refcache_set_chromosome(rc, chrm);
      } else if (strcmp(chrm, tsv_field(tsv, 0)) != 0) {
         fprintf(
            stderr,
            "[%s:%d] Error, rectangle cannot cross chromosomes.\n",
            __func__, __LINE__);
         exit(1);
      }

      // compute padding
      int p, pad = 0;
      for (p = region_beg; p < read_beg; ++pad)
         p = refcache_next_cg(rc, p) + 1;

      read_t *r = next_ref_read_v(reads);
      r->seq = (const kstring_t) {0};
      while(pad--) kputc('N', &r->seq);
      kputs(tsv->fields[5], &r->seq);
      r->other = (const kstring_t) {0};
      kputs(tsv->line, &r->other);

      // update region width
      if (region_width < 0 || (unsigned) region_width < r->seq.l)
         region_width = r->seq.l;
   }
   free(chrm);

   FILE *out = wzopen_out(out_fn);
   unsigned i;
   for (i = 0; i < reads->size; ++i) {
      read_t *r = ref_read_v(reads, i);
      while (r->seq.l < (unsigned) region_width) kputc('N', &r->seq);
      fputs(r->other.s, out);
      fputc('\t', out);
      fputs(r->seq.s, out);
      fputc('\n', out);
      free(r->seq.s);
      free(r->other.s);
   }

   tsv_close(tsv);
   free_read_v(reads);
   free_refcache(rc);

   return 0;
}

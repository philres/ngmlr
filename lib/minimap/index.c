#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "minimap.h"
#include "kvec.h"
#include "khash.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

mm_idx_t *mm_idx_init(int w, int k, int b)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b;
	mi->max_occ = UINT32_MAX;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	int i;
	if (mi == 0) return;
	for (i = 0; i < 1<<mi->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	free(mi->B);
	if (mi->name)
		for (i = 0; i < mi->n; ++i) free(mi->name[i]);
	free(mi->len); free(mi->name);
	free(mi);
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) {
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

uint32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return UINT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

void mm_idx_set_max_occ(mm_idx_t *mi, float f)
{
	mi->freq_thres = f;
	mi->max_occ = mm_idx_cal_max_occ(mi, f);
}

/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>mi->b<<1, &absent);
			assert(absent && j - start_a == n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == start_p);

	// deallocate and clear b->a
	free(b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	int tbatch_size, n_processed, keep_name;
	bseq_file_t *fp;
	uint64_t ibatch_size, n_read;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	bseq1_t *seq;
	mm128_v a;
} step_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x&mask].a;
		kv_push(mm128_t, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->n_read > p->ibatch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = bseq_read(p->fp, p->tbatch_size, &s->n_seq);
		if (s->seq) {
			uint32_t old_m = p->mi->n, m, n;
			assert((uint64_t)p->n_processed + s->n_seq <= INT32_MAX);
			m = n = p->mi->n + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m) {
				if (p->keep_name)
					p->mi->name = (char**)realloc(p->mi->name, m * sizeof(char*));
				p->mi->len = (int*)realloc(p->mi->len, m * sizeof(int));
			}
			for (i = 0; i < s->n_seq; ++i) {
				if (p->keep_name) {
					assert(strlen(s->seq[i].name) <= 254);
					p->mi->name[p->mi->n] = strdup(s->seq[i].name);
				}
				p->mi->len[p->mi->n++] = s->seq[i].l_seq;
				s->seq[i].rid = p->n_processed++;
				p->n_read += s->seq[i].l_seq;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			bseq1_t *t = &s->seq[i];
			mm_sketch(t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, &s->a);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		free(s->a.a); free(s);
	}
    return 0;
}

mm_idx_t *mm_idx_gen(bseq_file_t *fp, int w, int k, int b, int tbatch_size, int n_threads, uint64_t ibatch_size, int keep_name)
{
	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.tbatch_size = tbatch_size;
	pl.keep_name = keep_name;
	pl.ibatch_size = ibatch_size;
	pl.fp = fp;
	if (pl.fp == 0) return 0;
	pl.mi = mm_idx_init(w, k, b);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post(pl.mi, n_threads);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int n_threads) // a simpler interface
{
	bseq_file_t *fp;
	mm_idx_t *mi;
	fp = bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, MM_IDX_DEF_B, 1<<18, n_threads, UINT64_MAX, 1);
	mm_idx_set_max_occ(mi, 0.001);
	bseq_close(fp);
	return mi;
}

/*************
 * index I/O *
 *************/

#define MM_IDX_MAGIC "MMI\1"

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint32_t x[6];
	int i;
	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n, x[4] = mi->name? 1 : 0, x[5] = mi->max_occ;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 6, fp);
	fwrite(&mi->freq_thres, sizeof(float), 1, fp);
	fwrite(mi->len, 4, mi->n, fp);
	if (mi->name) {
		for (i = 0; i < mi->n; ++i) {
			uint8_t l;
			l = strlen(mi->name[i]);
			fwrite(&l, 1, 1, fp);
			fwrite(mi->name[i], 1, l, fp);
		}
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	int i;
	char magic[4];
	uint32_t x[6];
	mm_idx_t *mi;
	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 6, fp) != 6) return 0;
	mi = mm_idx_init(x[0], x[1], x[2]);
	mi->n = x[3], mi->max_occ = x[5];
	fread(&mi->freq_thres, sizeof(float), 1, fp);
	mi->len = (int32_t*)malloc(mi->n * 4);
	fread(mi->len, 4, mi->n, fp);
	if (x[4]) { // has names
		mi->name = (char**)calloc(mi->n, sizeof(char*));
		for (i = 0; i < mi->n; ++i) {
			uint8_t l;
			fread(&l, 1, 1, fp);
			mi->name[i] = (char*)malloc(l + 1);
			fread(mi->name[i], 1, l, fp);
			mi->name[i][l] = 0;
		}
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		fread(&size, 4, 1, fp);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, fp);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	return mi;
}

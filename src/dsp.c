#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "motif.h"

void *arena_alloc(Arena *a, size_t count, size_t size) {
    if (a->oom) return NULL;
    if (count && size > SIZE_MAX / count) { a->oom = 1; return NULL; }
    if (a->count == a->cap) {
        size_t cap = a->cap ? a->cap * 2 : 64;
        void **p = realloc(a->ptrs, cap * sizeof *p);
        if (!p) { a->oom = 1; return NULL; }
        a->ptrs = p;
        a->cap = cap;
    }
    void *mem = calloc(count ? count : 1, size ? size : 1);
    if (!mem) { a->oom = 1; return NULL; }
    a->ptrs[a->count++] = mem;
    return mem;
}

void *arena_grow(Arena *a, void *ptr, size_t count, size_t size) {
    if (!ptr) return arena_alloc(a, count, size);
    if (a->oom) return NULL;
    if (count && size > SIZE_MAX / count) { a->oom = 1; return NULL; }
    for (size_t i = a->count; i-- > 0;) {
        if (a->ptrs[i] != ptr) continue;
        void *mem = realloc(ptr, count * size);
        if (!mem) { a->oom = 1; return NULL; }
        a->ptrs[i] = mem;
        return mem;
    }
    abort(); /* not from this arena: a programming error */
}

size_t arena_mark(const Arena *a) { return a->count; }

void arena_release(Arena *a, size_t mark) {
    while (a->count > mark) free(a->ptrs[--a->count]);
}

void arena_free_all(Arena *a) {
    arena_release(a, 0);
    free(a->ptrs);
    a->ptrs = NULL;
    a->cap = 0;
}

void *vec_reserve(Arena *a, void *data, size_t *cap, size_t need, size_t elem) {
    if (need <= *cap) return data;
    size_t c = *cap ? *cap : 256;
    while (c < need) c *= 2;
    void *grown = arena_grow(a, data, c, elem);
    if (grown) *cap = c;
    return grown;
}

int rfft_init(Arena *a, Rfft *p, int n) {
    int m = n / 2;
    p->n = n;
    p->cos_t = arena_alloc(a, (size_t)m, sizeof(double));
    p->sin_t = arena_alloc(a, (size_t)m, sizeof(double));
    p->post_c = arena_alloc(a, (size_t)m + 1, sizeof(double));
    p->post_s = arena_alloc(a, (size_t)m + 1, sizeof(double));
    p->rev = arena_alloc(a, (size_t)m, sizeof(int));
    p->window = arena_alloc(a, (size_t)n, sizeof(float));
    p->zr = arena_alloc(a, (size_t)m, sizeof(double));
    p->zi = arena_alloc(a, (size_t)m, sizeof(double));
    if (a->oom) return -1;
    for (int j = 0; j < m; j++) {
        p->cos_t[j] = cos(2 * M_PI * j / m);
        p->sin_t[j] = sin(2 * M_PI * j / m);
    }
    for (int k = 0; k <= m; k++) {
        p->post_c[k] = cos(2 * M_PI * k / n);
        p->post_s[k] = sin(2 * M_PI * k / n);
    }
    p->post_c[m] = -1.0;
    p->post_s[m] = 0.0;
    int bits = 0;
    while ((1 << bits) < m) bits++;
    for (int k = 0; k < m; k++) {
        int r = 0;
        for (int b = 0; b < bits; b++) r |= ((k >> b) & 1) << (bits - 1 - b);
        p->rev[k] = r;
    }
    /* np.hanning(n + 1)[:-1], in NumPy's own formula. */
    p->window_sum = 0;
    for (int k = 0; k < n; k++) {
        double w = 0.5 + 0.5 * cos(M_PI * (2.0 * k - n) / n);
        p->window[k] = (float)w;
        p->window_sum += w;
    }
    return 0;
}

static void fft(Rfft *p) {
    int m = p->n / 2;
    double *re = p->zr, *im = p->zi;
    for (int size = 2; size <= m; size *= 2) {
        int half = size / 2, step = m / size;
        for (int start = 0; start < m; start += size) {
            for (int j = 0; j < half; j++) {
                double wr = p->cos_t[j * step], wi = -p->sin_t[j * step];
                int u = start + j, v = u + half;
                double tr = re[v] * wr - im[v] * wi, ti = re[v] * wi + im[v] * wr;
                re[v] = re[u] - tr;
                im[v] = im[u] - ti;
                re[u] += tr;
                im[u] += ti;
            }
        }
    }
}

void frame_power(Rfft *p, const float *x, size_t start, double *power) {
    int m = p->n / 2;
    const float *w = p->window;
    for (int k = 0; k < m; k++) {
        /* float32 product, as NumPy multiplies the float32 frame by the float32 window */
        p->zr[p->rev[k]] = (float)(x[start + 2 * k] * w[2 * k]);
        p->zi[p->rev[k]] = (float)(x[start + 2 * k + 1] * w[2 * k + 1]);
    }
    fft(p);
    /* Unpack the n/2-point complex FFT of (even, odd) samples into the real spectrum. */
    for (int k = 0; k <= m; k++) {
        double ar = p->zr[k % m], ai = p->zi[k % m];
        double br = p->zr[(m - k) % m], bi = p->zi[(m - k) % m];
        double er = (ar + br) / 2, ei = (ai - bi) / 2;
        double or_ = (ai + bi) / 2, oi = -(ar - br) / 2;
        double c = p->post_c[k], s = p->post_s[k];
        double xr = er + c * or_ + s * oi, xi = ei + c * oi - s * or_;
        power[k] = xr * xr + xi * xi;
    }
}

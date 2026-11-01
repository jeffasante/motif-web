/* Port of motif/matching.py: shared-audio segments and musical similarity. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "motif.h"

#define SEGMENT_GAP_SECONDS 2.0

static int by_hash(const void *pa, const void *pb) {
    const Hash *a = pa, *b = pb;
    if (a->hash != b->hash) return a->hash < b->hash ? -1 : 1;
    return (a->frame > b->frame) - (a->frame < b->frame);
}

int sort_hashes(Arena *a, const Hashes *in, Hashes *out) {
    out->items = arena_alloc(a, in->n, sizeof(Hash));
    if (!out->items) return -1;
    memcpy(out->items, in->items, in->n * sizeof(Hash));
    out->n = in->n;
    qsort(out->items, out->n, sizeof(Hash), by_hash);
    return 0;
}

/* First index in sorted[0..n) whose hash is >= h. */
static size_t lower_bound(const Hashes *sorted, uint32_t h) {
    size_t lo = 0, hi = sorted->n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (sorted->items[mid].hash < h) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

typedef struct {
    int32_t q, r;
    int64_t offset;
    size_t index;
} Pair;

static int by_q(const void *pa, const void *pb) {
    const Pair *a = pa, *b = pb;
    if (a->q != b->q) return a->q < b->q ? -1 : 1;
    return (a->index > b->index) - (a->index < b->index);
}

static int by_int64(const void *pa, const void *pb) {
    int64_t a = *(const int64_t *)pa, b = *(const int64_t *)pb;
    return (a > b) - (a < b);
}

typedef struct {
    int64_t value;
    size_t window;
} Candidate;

static int by_window(const void *pa, const void *pb) {
    const Candidate *a = pa, *b = pb;
    if (a->window != b->window) return a->window > b->window ? -1 : 1;
    return (a->value > b->value) - (a->value < b->value); /* stable argsort of unique (sorted) values */
}

static int by_strength(const void *pa, const void *pb) {
    const Segment *a = pa, *b = pb;
    if (a->aligned != b->aligned) return a->aligned > b->aligned ? -1 : 1;
    return (a->query_start > b->query_start) - (a->query_start < b->query_start);
}

/* find_segments: group pairs by offset (+-tol), strongest first, split into time runs. */
static int find_segments(Arena *a, Pair *pairs, size_t n, double frame_seconds, const MatchingConfig *c,
                         Segment *segments, int *n_segments) {
    *n_segments = 0;
    if (n == 0) return 0;
    size_t mark = arena_mark(a);
    int64_t *sorted = arena_alloc(a, n, sizeof *sorted);
    Candidate *cand = arena_alloc(a, n, sizeof *cand);
    size_t *counts = arena_alloc(a, n, sizeof *counts);
    unsigned char *used = arena_alloc(a, n, 1);
    Pair *members = arena_alloc(a, n, sizeof *members);
    if (!sorted || !cand || !counts || !used || !members) return -1;
    for (size_t i = 0; i < n; i++) sorted[i] = pairs[i].offset;
    qsort(sorted, n, sizeof *sorted, by_int64);
    size_t u = 0;
    for (size_t i = 0; i < n; i++) {
        if (u && cand[u - 1].value == sorted[i]) { counts[u - 1]++; continue; }
        cand[u].value = sorted[i];
        counts[u++] = 1;
    }
    int64_t tol = c->offset_tolerance_frames;
    size_t lo = 0, hi = 0, sum = 0;
    for (size_t k = 0; k < u; k++) {
        while (hi < u && cand[hi].value <= cand[k].value + tol) sum += counts[hi++];
        while (cand[lo].value < cand[k].value - tol) sum -= counts[lo++];
        cand[k].window = sum;
    }
    qsort(cand, u, sizeof *cand, by_window);
    long gap = lrint(SEGMENT_GAP_SECONDS / frame_seconds);
    if (gap < 1) gap = 1;
    int count = 0;
    Segment found[64];
    for (size_t k = 0; k < u && count < c->max_segments; k++) {
        if (cand[k].window < (size_t)c->min_aligned_hashes) break;
        int64_t centre = cand[k].value;
        size_t m = 0;
        for (size_t i = 0; i < n; i++) {
            int64_t d = pairs[i].offset - centre;
            if (!used[i] && (d < 0 ? -d : d) <= tol) {
                members[m] = pairs[i];
                members[m++].index = i;
            }
        }
        if (m < (size_t)c->min_aligned_hashes) continue;
        qsort(members, m, sizeof *members, by_q);
        for (size_t first = 0; first < m;) {
            size_t last = first + 1;
            while (last < m && members[last].q - members[last - 1].q <= gap) last++;
            if (last - first >= (size_t)c->min_aligned_hashes && count < 64) {
                int32_t q0 = members[first].q, q1 = members[last - 1].q;
                int32_t r0 = members[first].r, r1 = members[first].r;
                for (size_t i = first; i < last; i++) {
                    used[members[i].index] = 1;
                    if (members[i].r < r0) r0 = members[i].r;
                    if (members[i].r > r1) r1 = members[i].r;
                }
                found[count++] = (Segment){q0 * frame_seconds, q1 * frame_seconds, r0 * frame_seconds,
                                           r1 * frame_seconds, (double)centre * frame_seconds, (int)(last - first)};
            }
            first = last;
        }
    }
    qsort(found, (size_t)count, sizeof *found, by_strength);
    if (count > c->max_segments) count = c->max_segments;
    if (count > MAX_SEGMENTS) count = MAX_SEGMENTS;
    memcpy(segments, found, (size_t)count * sizeof *found);
    *n_segments = count;
    arena_release(a, mark);
    return 0;
}

static int shared_audio(Arena *a, const Hashes *query, const Hashes *ref, const FingerprintConfig *fp,
                        const MatchingConfig *c, double ratio, SharedAudio *out) {
    size_t mark = arena_mark(a);
    memset(out, 0, sizeof *out);
    out->speed_ratio = ratio;
    out->query_hashes = query->n;
    size_t n_pairs = 0;
    for (size_t i = 0; i < query->n; i++) {
        size_t j = lower_bound(ref, query->items[i].hash), first = j;
        while (j < ref->n && ref->items[j].hash == query->items[i].hash) j++;
        if (j > first) out->matching_hashes++;
        n_pairs += j - first;
    }
    Pair *pairs = arena_alloc(a, n_pairs, sizeof *pairs);
    if (!pairs) return -1;
    size_t p = 0;
    for (size_t i = 0; i < query->n; i++) {
        for (size_t j = lower_bound(ref, query->items[i].hash); j < ref->n && ref->items[j].hash == query->items[i].hash; j++) {
            pairs[p].q = query->items[i].frame;
            pairs[p].r = ref->items[j].frame;
            pairs[p].offset = (int64_t)pairs[p].r - pairs[p].q;
            p++;
        }
    }
    double frame_seconds = (double)fp->hop_length / SAMPLE_RATE;
    if (find_segments(a, pairs, n_pairs, frame_seconds, c, out->segments, &out->n_segments)) return -1;
    for (int i = 0; i < out->n_segments; i++) {
        /* Query frames were scaled into reference time; report them in the query's own time. */
        out->segments[i].query_start /= ratio;
        out->segments[i].query_end /= ratio;
    }
    out->detected = out->n_segments > 0;
    arena_release(a, mark);
    return 0;
}

static int strongest(const SharedAudio *s) { return s->n_segments ? s->segments[0].aligned : 0; }

int shared_audio_any_speed(Arena *a, const Peaks *query, const Hashes *reference_sorted,
                           const FingerprintConfig *fp, const MatchingConfig *c, SharedAudio *out) {
    SharedAudio result;
    for (int i = 0; i < c->n_speed_ratios; i++) {
        size_t mark = arena_mark(a);
        Hashes q;
        if (hash_peaks(a, query, fp, c->speed_ratios[i], &q)) return -1;
        if (shared_audio(a, &q, reference_sorted, fp, c, c->speed_ratios[i], &result)) return -1;
        if (i == 0 || strongest(&result) > strongest(out)) *out = result;
        arena_release(a, mark);
    }
    return 0;
}

int align_hashes(Arena *a, const Hashes *query, const Hashes *ref, int tolerance, Alignment *out) {
    memset(out, 0, sizeof *out);
    out->matched = arena_alloc(a, query->n + 1, 1);
    if (!out->matched) return -1;
    for (size_t i = 0; i < query->n; i++) {
        size_t j = lower_bound(ref, query->items[i].hash), first = j;
        while (j < ref->n && ref->items[j].hash == query->items[i].hash) j++;
        out->matched[i] = j > first;
        out->matching += j > first;
        out->n_pairs += j - first;
    }
    out->q = arena_alloc(a, out->n_pairs + 1, sizeof *out->q);
    out->offset = arena_alloc(a, out->n_pairs + 1, sizeof *out->offset);
    int64_t *sorted = arena_alloc(a, out->n_pairs + 1, sizeof *sorted);
    if (!out->q || !out->offset || !sorted) return -1;
    size_t p = 0;
    for (size_t i = 0; i < query->n; i++) {
        for (size_t j = lower_bound(ref, query->items[i].hash); j < ref->n && ref->items[j].hash == query->items[i].hash; j++) {
            out->q[p] = query->items[i].frame;
            out->offset[p] = (int64_t)ref->items[j].frame - query->items[i].frame;
            sorted[p] = out->offset[p];
            p++;
        }
    }
    qsort(sorted, p, sizeof *sorted, by_int64);
    size_t lo = 0;
    for (size_t hi = 0; hi < p; hi++) {
        while (sorted[hi] - sorted[lo] > 2 * (int64_t)tolerance) lo++;
        if (hi - lo + 1 > out->best_count) {
            out->best_count = hi - lo + 1;
            out->best_offset = sorted[lo] + (sorted[hi] - sorted[lo]) / 2;
        }
    }
    return 0;
}

static int signed_semitones(int value) {
    int v = ((value % 12) + 12) % 12;
    return v > 6 ? v - 12 : v;
}

/* Per frame: subtract the mean over pitch classes and scale to unit norm (zero if flat). */
static double *standardize(Arena *a, const Chroma *ch) {
    int n = ch->frames;
    double *out = arena_alloc(a, (size_t)12 * n, sizeof(double));
    if (!out) return NULL;
    for (int t = 0; t < n; t++) {
        double mean = 0, norm = 0;
        for (int k = 0; k < 12; k++) mean += ch->data[k * n + t];
        mean /= 12;
        for (int k = 0; k < 12; k++) {
            double v = ch->data[k * n + t] - mean;
            out[k * n + t] = v;
            norm += v * v;
        }
        norm = sqrt(norm);
        for (int k = 0; k < 12; k++) out[k * n + t] = norm > 0 ? out[k * n + t] / norm : 0;
    }
    return out;
}

/* Resample the time axis by ratio (linear interpolation); returns frame-major [t*12+k]. */
static double *stretch(Arena *a, const double *x, int n, double ratio, int *n_out) {
    long m = lrint(n * ratio);
    if (m < 1) m = 1;
    double *out = arena_alloc(a, (size_t)m * 12, sizeof(double));
    if (!out) return NULL;
    for (long t = 0; t < m; t++) {
        double src = (double)t / ratio;
        long i0 = (long)floor(src);
        if (i0 < 0) i0 = 0;
        if (i0 > n - 1) i0 = n - 1;
        long i1 = i0 + 1 > n - 1 ? n - 1 : i0 + 1;
        double w = src - (double)i0;
        for (int k = 0; k < 12; k++) out[t * 12 + k] = (1 - w) * x[k * n + i0] + w * x[k * n + i1];
    }
    *n_out = (int)m;
    return out;
}

int musical_similarity(Arena *a, const Chroma *query, const Chroma *reference, const MusicalConfig *mc,
                       const MatchingConfig *c, Musical *out) {
    memset(out, 0, sizeof *out);
    for (int k = 0; k < 12; k++) out->shift_best[k] = -INFINITY;
    for (int k = 0; k < MAX_TEMPO_RATIOS; k++) out->ratio_best[k] = -INFINITY;
    double frame = mc->chroma_frame_seconds;
    int width = (int)lrint(c->melodic_window_seconds / frame);
    int hop = (int)lrint(c->melodic_window_hop_seconds / frame);
    if (hop < 1) hop = 1;
    int m = reference->frames;
    if (m < width) {
        out->short_reference = 1;
        return 0;
    }
    size_t mark = arena_mark(a);
    double *ref = standardize(a, reference), *qs = standardize(a, query);
    double *refT = arena_alloc(a, (size_t)m * 12, sizeof(double));
    if (!ref || !qs || !refT) return -1;
    for (int j = 0; j < m; j++)
        for (int k = 0; k < 12; k++) refT[j * 12 + k] = ref[k * m + j];
    int cols = m - width + 1;
    int ring = width / hop + 2;
    double best_evidence = 0;
    int have_best = 0;
    for (int ri = 0; ri < c->n_melodic_tempo_ratios; ri++) {
        double ratio = c->melodic_tempo_ratios[ri];
        size_t ratio_mark = arena_mark(a);
        int n = 0;
        double *q = stretch(a, qs, query->frames, ratio, &n);
        if (!q) return -1;
        if (n < width) {
            arena_release(a, ratio_mark);
            continue;
        }
        int rows_total = n - width + 1, n_starts = (rows_total + hop - 1) / hop;
        int last = (n_starts - 1) * hop, n_offsets = last + cols;
        /* C_k[l]: cumulative diagonal sums of sim, exactly as _diagonal_window_means builds them. */
        double *prev = arena_alloc(a, (size_t)m + 1, sizeof(double));
        double *cur = arena_alloc(a, (size_t)m + 1, sizeof(double));
        double *snap = arena_alloc(a, (size_t)ring * (m + 1), sizeof(double));
        double *row = arena_alloc(a, (size_t)cols, sizeof(double));
        double *tol = arena_alloc(a, (size_t)cols, sizeof(double));
        int *run = arena_alloc(a, (size_t)n_offsets, sizeof(int));
        double *total = arena_alloc(a, (size_t)n_offsets, sizeof(double));
        if (!prev || !cur || !snap || !row || !tol || !run || !total) return -1;
        for (int shift = 0; shift < 12; shift++) {
            memset(prev, 0, (size_t)(m + 1) * sizeof(double));
            memset(run, 0, (size_t)n_offsets * sizeof(int));
            memset(total, 0, (size_t)n_offsets * sizeof(double));
            memcpy(snap, prev, (size_t)(m + 1) * sizeof(double)); /* C_0 for start 0 */
            for (int k = 1; k <= n; k++) {
                /* row k-1 of sim = roll(q, shift).T @ ref */
                const double *qv = q + (size_t)(k - 1) * 12;
                double rolled[12];
                for (int p = 0; p < 12; p++) rolled[p] = qv[(p - shift + 12) % 12];
                cur[0] = 0;
                for (int l = 1; l <= m; l++) {
                    const double *rv = refT + (size_t)(l - 1) * 12;
                    double s = 0;
                    for (int p = 0; p < 12; p++) s += rolled[p] * rv[p];
                    cur[l] = s + prev[l - 1];
                }
                double *swap = prev;
                prev = cur;
                cur = swap; /* prev now holds C_k */
                if (k % hop == 0 && k <= last)
                    memcpy(snap + (size_t)((k / hop) % ring) * (m + 1), prev, (size_t)(m + 1) * sizeof(double));
                int s = k - width;
                if (s < 0 || s % hop != 0 || s > last) continue;
                int w = s / hop;
                const double *base = snap + (size_t)(w % ring) * (m + 1);
                double *sb = &out->shift_best[(12 - shift) % 12];
                double *rb = &out->ratio_best[ri < MAX_TEMPO_RATIOS ? ri : MAX_TEMPO_RATIOS - 1];
                for (int j = 0; j < cols; j++) {
                    row[j] = (prev[j + width] - base[j]) / width;
                    if (row[j] > *sb) *sb = row[j];
                    if (row[j] > *rb) *rb = row[j];
                }
                for (int j = 0; j < cols; j++) {
                    double v = row[j];
                    if (j > 0 && row[j - 1] > v) v = row[j - 1];
                    if (j + 1 < cols && row[j + 1] > v) v = row[j + 1];
                    tol[j] = v;
                }
                int any = 0;
                for (int o = 0; o < n_offsets; o++) {
                    int col = s + o - last;
                    double score = col >= 0 && col < cols ? tol[col] : -INFINITY;
                    if (score >= c->melodic_min_similarity) {
                        run[o]++;
                        total[o] += score;
                        if (run[o] >= c->melodic_min_run_windows) any = 1;
                    } else {
                        run[o] = 0;
                        total[o] = 0;
                    }
                }
                if (!any) continue;
                /* Rank by matched query seconds x mean similarity. */
                int bk = -1;
                double be = 0, bm = 0;
                for (int o = 0; o < n_offsets; o++) {
                    double mean = 0, ev = 0;
                    if (run[o] >= c->melodic_min_run_windows) {
                        mean = total[o] / run[o];
                        ev = ((run[o] - 1) * hop + width) / ratio * frame * mean;
                    }
                    if (bk < 0 || ev > be) { bk = o; be = ev; bm = mean; }
                }
                if (!have_best || be > best_evidence) {
                    have_best = 1;
                    best_evidence = be;
                    int length = run[bk], d = bk - last;
                    int q0 = (w - length + 1) * hop, q1 = s + width;
                    out->detected = 1;
                    out->windows = length;
                    out->mean_similarity = bm;
                    out->transposition = signed_semitones(-shift);
                    out->tempo_ratio = ratio;
                    out->query_start = q0 / ratio * frame;
                    out->query_end = q1 / ratio * frame;
                    out->reference_start = (q0 + d) * frame;
                    out->reference_end = (q1 + d) * frame;
                }
            }
        }
        arena_release(a, ratio_mark);
    }
    arena_release(a, mark);
    return 0;
}

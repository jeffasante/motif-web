/* Motif C engine: a port of motif.features and motif.matching for `compare`,
 * with the Python defaults built in. Internal interface; the export is motif_compare. */
#ifndef MOTIF_H
#define MOTIF_H

#include <stddef.h>
#include <stdint.h>

#define MOTIF_VERSION "0.1.0"
#define SAMPLE_RATE 22050
#define MAX_SEGMENTS 5
#define MAX_TEMPO_RATIOS 8

/* Every allocation of one motif_compare call, freed together. On failure `oom` is set
 * and NULL returned; callers return -1 and the entry point reports out of memory. */
typedef struct {
    void **ptrs;
    size_t count, cap;
    int oom;
} Arena;

void *arena_alloc(Arena *a, size_t count, size_t size); /* zeroed */
void *arena_grow(Arena *a, void *ptr, size_t count, size_t size);
size_t arena_mark(const Arena *a);
void arena_release(Arena *a, size_t mark); /* frees everything allocated after mark */
void arena_free_all(Arena *a);
/* Returns data, moved if it had to grow, with room for `need` elements; NULL on failure. */
void *vec_reserve(Arena *a, void *data, size_t *cap, size_t need, size_t elem);

typedef struct {
    int n;
    double *cos_t, *sin_t;   /* complex FFT of n/2 points */
    double *post_c, *post_s; /* real-spectrum unpacking, k = 0..n/2 */
    int *rev;
    float *window;           /* periodic Hann, rounded to float32 as NumPy does */
    double window_sum;       /* of the float64 window */
    double *zr, *zi;
} Rfft;

int rfft_init(Arena *a, Rfft *p, int n);
/* power[0..n/2] = |rfft(x[start:start+n] * window)|^2 */
void frame_power(Rfft *p, const float *x, size_t start, double *power);

static inline int frame_count(size_t n, int n_fft, int hop) {
    return n < (size_t)n_fft ? 0 : 1 + (int)((n - (size_t)n_fft) / (size_t)hop);
}

typedef struct {
    int n_fft, hop_length;
    double f_min, f_max, min_peak_dbfs;
    int neighborhood_frames, neighborhood_bins, peaks_per_second, fan_out;
    int target_min_frames, target_max_frames, target_max_bins;
} FingerprintConfig;

typedef struct {
    int chroma_n_fft, chroma_hop_length;
    double chroma_f_min, chroma_f_max, chroma_frame_seconds, key_min_correlation;
    double tempo_min_bpm, tempo_max_bpm, tempo_min_strength;
} MusicalConfig;

typedef struct {
    int min_aligned_hashes, offset_tolerance_frames, max_segments;
    const double *speed_ratios;
    int n_speed_ratios;
    double melodic_window_seconds, melodic_window_hop_seconds, melodic_min_similarity;
    int melodic_min_run_windows;
    const double *melodic_tempo_ratios;
    int n_melodic_tempo_ratios;
} MatchingConfig;

typedef struct {
    int32_t frame, bin, bucket;
    float level;
} PeakRow;

typedef struct {
    PeakRow *rows; /* by frame, then bin */
    size_t n;
    int *found_per_bucket; /* peaks found in each bucket before the peaks_per_second cap */
    int n_buckets, frames_per_bucket;
} Peaks;

typedef struct {
    uint32_t hash;
    int32_t frame;
} Hash;

typedef struct {
    Hash *items;
    size_t n;
} Hashes;

typedef struct {
    float *data; /* data[pc * frames + t] */
    int frames;
} Chroma;

typedef struct {
    int has_scores, known, best, second; /* best/second: 0-11 major, 12-23 minor */
    double best_score, second_score;
} KeyEstimate;

enum { TEMPO_NO_ONSETS, TEMPO_WEAK, TEMPO_KNOWN };

typedef struct {
    int state, n_alternatives;
    double bpm, strength, alternatives[2];
} TempoEstimate;

typedef struct {
    double query_start, query_end, reference_start, reference_end, offset;
    int aligned;
} Segment;

typedef struct {
    int detected, n_segments;
    double speed_ratio;
    size_t matching_hashes, query_hashes;
    Segment segments[MAX_SEGMENTS];
} SharedAudio;

typedef struct {
    int detected, short_reference, windows, transposition;
    double mean_similarity, tempo_ratio, query_start, query_end, reference_start, reference_end;
    /* Best single-window similarity per transposition (index = semitones mod 12) and per
     * tempo ratio; -INFINITY where nothing was scored. */
    double shift_best[12], ratio_best[MAX_TEMPO_RATIOS];
} Musical;

/* Every (query hash, reference hash) pair with equal hash values, for the diagnostics. */
typedef struct {
    size_t n_pairs, matching;  /* matching: query hashes with at least one reference match */
    int32_t *q;                /* query frame of each pair */
    int64_t *offset;           /* reference frame - query frame */
    unsigned char *matched;    /* per query hash */
    int64_t best_offset;
    size_t best_count;         /* pairs within +-tolerance of best_offset */
} Alignment;

int align_hashes(Arena *a, const Hashes *query, const Hashes *reference_sorted, int tolerance, Alignment *out);

extern const char *const PITCH_CLASSES[12];

int find_peaks(Arena *a, const float *x, size_t n, const FingerprintConfig *c, Peaks *out);
int hash_peaks(Arena *a, const Peaks *peaks, const FingerprintConfig *c, double speed_ratio, Hashes *out);
/* As hash_peaks, with frequencies divided by freq_ratio and times multiplied by time_ratio. */
int hash_peaks_scaled(Arena *a, const Peaks *peaks, const FingerprintConfig *c, double freq_ratio, double time_ratio,
                      Hashes *out);
int chroma(Arena *a, const float *x, size_t n, const MusicalConfig *c, Chroma *out);
void estimate_key(const Chroma *ch, const MusicalConfig *c, KeyEstimate *out);
int estimate_tempo(Arena *a, const float *x, size_t n, const MusicalConfig *c, TempoEstimate *out);

int sort_hashes(Arena *a, const Hashes *in, Hashes *out);
int shared_audio_any_speed(Arena *a, const Peaks *query, const Hashes *reference_sorted,
                           const FingerprintConfig *fp, const MatchingConfig *c, SharedAudio *out);
int musical_similarity(Arena *a, const Chroma *query, const Chroma *reference, const MusicalConfig *mc,
                       const MatchingConfig *c, Musical *out);

#endif

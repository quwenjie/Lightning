#pragma once
// brakedown.hpp — Brakedown linear-time code (mcl BN254 Fr)
// Parameters: α=0.1195, r=1.42, c=6, d=33 (Figure 2 row 1)

#include <mcl/bn.hpp>
#include <vector>
#include <random>
#include <cassert>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>

using F = mcl::bn::Fr;

// ==================== NTT utilities ====================
F pow_F(F base, uint64_t exp);
F getRoot(int n);
void ntt(std::vector<F>& a, const F& root);
void intt(std::vector<F>& a, const F& root);
void rs_ntt(const F* msg, int k, F* out, int n);

// ==================== Sparse Matrix ====================
struct SparseMatrix {
    int rows, cols_count, sp;
    std::vector<int> col_idx;   // flat: rows * sp
    std::vector<F>   vals;      // flat: rows * sp

    void generate(int _rows, int _cols, int _sp, std::mt19937_64& rng);
    size_t mem_bytes() const;
};

void sparse_mul(const F* x, const SparseMatrix& M, F* y);

// ==================== RS base code ====================
struct RSCode {
    int k;      // message length
    int n;      // codeword length (power of 2)
    F root;     // primitive n-th root of unity

    void setup(int msg_len, double target_rate);
    void encode(const F* x, F* out) const;
};

// ==================== Brakedown Code ====================
// Enc(x) = [x | z | v]
//   y = x·A   (n → αn)
//   z = Enc(y) (recursive, αn → αrn)
//   v = z·B   (αrn → (r-1-rα)n)

struct BrakedownLevel {
    int n, an, arn, tail;
    SparseMatrix A, B;
};

struct BrakedownCode {
    static constexpr double ALPHA = 0.1195;
    static constexpr double R_PAR = 1.42;
    static constexpr int    C_SP  = 6;
    static constexpr int    D_SP  = 33;
    static constexpr int    BASE_N = (1 << 18);

    int n_;                            // top-level message length
    std::vector<BrakedownLevel> levels_;
    RSCode rs_base_;
    int cw_len_;                       // total codeword length

    // Setup for encoding messages of length msg_len.
    void setup(int msg_len, std::mt19937_64& rng);

    // Encode message x of length msg_len into out.
    // Returns actual codeword length written.
    int encode(const F* x, F* out) const;

    // Expected codeword length (computed at setup).
    int codeword_len() const { return cw_len_; }

private:
    int encode_recursive(const F* x, int msg_len, int depth, F* out) const;
    int compute_cw_len(int msg_len, int depth) const;
};

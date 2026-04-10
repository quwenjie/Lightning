// brakedown.cpp — Brakedown linear-time code implementation
#include "brakedown.hpp"
#include <iostream>

// ==================== NTT ====================
F pow_F(F base, uint64_t exp) {
    F r; r = 1;
    while (exp) {
        if (exp & 1) r *= base;
        base *= base;
        exp >>= 1;
    }
    return r;
}

F getRoot(int n) {
    F res;
    res = -F(1);
    if (n == 1) { res = 1; return res; }
    while (n > 1) {
        F::squareRoot(res, res);
        n >>= 1;
    }
    return res;
}

void ntt(std::vector<F>& a, const F& root) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        F wlen = pow_F(root, n / len);
        for (int i = 0; i < n; i += len) {
            F w; w = 1;
            for (int j = 0; j < len / 2; j++) {
                F u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j]           = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

void intt(std::vector<F>& a, const F& root) {
    F inv_root;
    F::inv(inv_root, root);
    ntt(a, inv_root);
    F inv_n;
    F::inv(inv_n, F(a.size()));
    for (auto& x : a) x *= inv_n;
}

void rs_ntt(const F* msg, int k, F* out, int n) {
    int N = 1;
    while (N < n) N <<= 1;
    std::vector<F> a(N);
    for (int i = 0; i < k; i++) a[i] = msg[i];
    for (int i = k; i < N; i++) a[i] = F(0);
    F r = getRoot(N);
    ntt(a, r);
    for (int i = 0; i < n; i++) out[i] = a[i];
}

// ==================== Sparse Matrix ====================
void SparseMatrix::generate(int _rows, int _cols, int _sp, std::mt19937_64& rng) {
    rows = _rows;
    cols_count = _cols;
    sp = std::min(_sp, _cols);
    size_t total = (size_t)rows * sp;
    col_idx.resize(total);
    vals.resize(total);

    std::vector<int> pool(cols_count);
    std::iota(pool.begin(), pool.end(), 0);

    for (int i = 0; i < rows; i++) {
        // partial Fisher-Yates
        for (int j = 0; j < sp; j++) {
            std::uniform_int_distribution<int> dist(j, cols_count - 1);
            int r = dist(rng);
            std::swap(pool[j], pool[r]);
        }
        size_t base = (size_t)i * sp;
        for (int j = 0; j < sp; j++) {
            col_idx[base + j] = pool[j];
            vals[base + j].setByCSPRNG();
            while (vals[base + j].isZero())
                vals[base + j].setByCSPRNG();
        }
    }
}

size_t SparseMatrix::mem_bytes() const {
    return col_idx.size() * sizeof(int) + vals.size() * sizeof(F);
}

void sparse_mul(const F* x, const SparseMatrix& M, F* y) {
    for (int j = 0; j < M.cols_count; j++) y[j] = 0;
    for (int i = 0; i < M.rows; i++) {
        if (x[i].isZero()) continue;
        size_t base = (size_t)i * M.sp;
        for (int j = 0; j < M.sp; j++) {
            y[M.col_idx[base + j]] += x[i] * M.vals[base + j];
        }
    }
}

// ==================== RS Code ====================
void RSCode::setup(int msg_len, double target_rate) {
    k = msg_len;
    int min_n = (int)std::ceil((double)k / target_rate);
    n = 1;
    while (n < min_n) n <<= 1;
    root = getRoot(n);
}

void RSCode::encode(const F* x, F* out) const {
    std::vector<F> poly(n);
    for (int i = 0; i < k; i++) poly[i] = x[i];
    for (int i = k; i < n; i++) poly[i] = F(0);
    ntt(poly, root);
    std::memcpy(out, poly.data(), n * sizeof(F));
}

// ==================== Brakedown Code ====================
void BrakedownCode::setup(int msg_len, std::mt19937_64& rng) {
    n_ = msg_len;
    levels_.clear();

    int cur = n_;
    while (cur >= BASE_N) {
        BrakedownLevel lv;
        lv.n    = cur;
        lv.an   = std::max(1, (int)std::round(ALPHA * cur));
        lv.arn  = std::max(1, (int)std::round(ALPHA * R_PAR * cur));
        lv.tail = std::max(1, (int)std::round((R_PAR - 1.0 - R_PAR * ALPHA) * cur));
        lv.A.generate(lv.n,   lv.an,  C_SP, rng);
        lv.B.generate(lv.arn, lv.tail, D_SP, rng);
        cur = lv.an;
        levels_.push_back(std::move(lv));
    }

    rs_base_.setup(cur, 1.0 / R_PAR);
    cw_len_ = compute_cw_len(n_, 0);
}

int BrakedownCode::compute_cw_len(int msg_len, int depth) const {
    if (msg_len < BASE_N)
        return rs_base_.n;
    const auto& lv = levels_[depth];
    int z_len = compute_cw_len(lv.an, depth + 1);
    if (z_len < lv.arn) z_len = lv.arn;
    return lv.n + z_len + lv.tail;
}

int BrakedownCode::encode(const F* x, F* out) const {
    return encode_recursive(x, n_, 0, out);
}

int BrakedownCode::encode_recursive(const F* x, int msg_len, int depth, F* out) const {
    if (msg_len < BASE_N) {
        rs_base_.encode(x, out);
        return rs_base_.n;
    }
    const auto& lv = levels_[depth];
    assert(msg_len == lv.n);

    // 1) y = x · A
    std::vector<F> y(lv.an);
    sparse_mul(x, lv.A, y.data());

    // 2) z = Enc(y) placed at out + n
    F* z_ptr = out + lv.n;
    int z_len = encode_recursive(y.data(), lv.an, depth + 1, z_ptr);
    if (z_len < lv.arn) {
        for (int i = z_len; i < lv.arn; i++) z_ptr[i] = F(0);
        z_len = lv.arn;
    }

    // 3) v = z · B
    F* v_ptr = out + lv.n + z_len;
    sparse_mul(z_ptr, lv.B, v_ptr);

    // 4) systematic: copy x to front
    std::memcpy(out, x, lv.n * sizeof(F));
    return lv.n + z_len + lv.tail;
}

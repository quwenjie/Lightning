// lightning.cpp — Lightning PCS with Brakedown inner code
//
// Compile:
//   g++ -O2 -std=c++17 -o lightning brakedown.cpp lightning.cpp -lmcl -lgmp -lcrypto
//
// Run:
//   ./lightning 24 0.004

#include "brakedown.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <numeric>
#include <openssl/sha.h>
#include <cassert>
#include <map>
#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <string>
#include <stdexcept>

using namespace mcl::bn;
using namespace std;
using namespace std::chrono;

// ==================== Lightning parameters ====================
struct Params {
    double alpha;
    int degree;
};

map<double, Params> table = {
    {0.003, {0.198, 10}},
    {0.004, {0.274, 9}},
    {0.005, {0.343, 9}},
    {0.006, {0.417,  9}},
    {0.007, {0.495,  9}},
    {0.008, {0.579,  9}},
};

double compute_t(double delta) {
    double x = 1.0 - delta / 3.0;
    double denom = -log2(x);
    return 90.0 / denom;
}

int get_best_col_num(int N, int t) {
    int best = 1;
    long long best_ps = (1LL << 62);
    for (int col = 1; col <= N - 1; col++) {
        int row = N - col;
        long long size = (1LL << row) * 1LL * t + 3LL * (1LL << col);
        if (size < best_ps) {
            best_ps = size;
            best = col;
        }
    }
    return best;
}

// ==================== Hash ====================
using Hash = array<uint8_t, 32>;

inline Hash sha256(const uint8_t* data, size_t len) {
    Hash out;
    SHA256(data, len, out.data());
    return out;
}

vector<int> sample_positions(int n, int t) {
    vector<int> all(n);
    iota(all.begin(), all.end(), 0);
    random_device rd;
    mt19937_64 rng(rd());
    vector<int> out;
    out.reserve(t);
    sample(all.begin(), all.end(), back_inserter(out), t, rng);
    return out;
}

// ==================== LDPC matrix (binary, D-row-regular) ====================
struct Random_matrix {
    int left_size, right_size, degree, row_number;
    vector<vector<int>> edges;

    Random_matrix(int lsize, int rsize, int deg, int rows)
        : left_size(lsize), right_size(rsize), degree(deg), row_number(rows)
    {
        edges.assign(left_size, vector<int>());
        vector<int> pool(right_size);
        iota(pool.begin(), pool.end(), 0);

        random_device rd;
        mt19937_64 rng(rd());

        for (int u = 0; u < left_size; u++) {
            for (int i = 0; i < degree; i++) {
                uniform_int_distribution<int> dist(i, right_size - 1);
                int j = dist(rng);
                swap(pool[i], pool[j]);
            }
            edges[u].assign(pool.begin(), pool.begin() + degree);
        }
    }

    void aggregate_Fr(const F* in, F* out) const {
        for (int u = 0; u < left_size; u++) {
            const auto& R = edges[u];
            for (int v : R) {
                out[v] += in[u];
            }
        }
    }

    void aggregate_batch_left_major(const uint64_t* in, F* out) const {
        vector<F> acc((size_t)right_size * row_number, F(0));

        for (int u = 0; u < left_size; u++) {
            const uint64_t* in_u = in + (size_t)u * row_number;
            const auto& R = edges[u];
            for (int v : R) {
                for (int i = 0; i < row_number; i++) {
                    acc[(size_t)v * row_number + i] += F(in_u[i]);
                }
            }
        }

        for (int v = 0; v < right_size; v++) {
            F* out_v = out + (size_t)v * row_number;
            for (int i = 0; i < row_number; i++) {
                out_v[i] = acc[(size_t)v * row_number + i];
            }
        }
    }
};

static void print_usage(const char* prog) {
    cerr << "Usage: " << prog << " <poly_var> <delta>\n";
    cerr << "Example:\n";
    cerr << "  " << prog << " 24 0.004\n";
    cerr << "Supported deltas:";
    for (const auto& kv : table) cerr << " " << kv.first;
    cerr << "\n";
}

int main(int argc, char** argv)
{
    if (argc != 3) {
        print_usage(argv[0]);
        return 1;
    }

    int poly_var;
    double delta;
    try {
        poly_var = stoi(argv[1]);
        delta = stod(argv[2]);
    } catch (const exception&) {
        cerr << "Invalid arguments. Expect: integer poly_var and floating-point delta.\n";
        print_usage(argv[0]);
        return 1;
    }

    if (poly_var <= 0 || poly_var >= 63) {
        cerr << "poly_var must be in 1..62.\n";
        return 1;
    }

    auto it = table.find(delta);
    if (it == table.end()) {
        cerr << "Unsupported delta: " << delta << "\n";
        print_usage(argv[0]);
        return 1;
    }

    const double alpha = it->second.alpha;
    const int degree   = it->second.degree;
    const int t        = static_cast<int>(compute_t(delta));

    const int best_col_log = get_best_col_num(poly_var, t);
    const int COLUMN_NUMBER = 1LL << best_col_log;
    const int ROW_NUMBER    = (1LL << poly_var) / COLUMN_NUMBER;
    const int right_size    = (int)ceil(alpha * COLUMN_NUMBER);

    initPairing(mcl::BN_SNARK1);

    mt19937_64 rng_bd(54321);
    BrakedownCode inner_code;
    inner_code.setup(right_size, rng_bd);
    const int second_len   = inner_code.codeword_len();
    const int codeword_len = COLUMN_NUMBER + second_len;

    cout << "=== Lightning PCS (Brakedown inner code) ===" << endl;
    cout << "poly_var=" << poly_var << "  delta=" << delta << endl;
    cout << "alpha=" << alpha << "  degree=" << degree << "  t=" << t << endl;
    cout << "COLUMN_NUMBER=" << COLUMN_NUMBER << "  ROW_NUMBER=" << ROW_NUMBER << endl;
    cout << "right_size=" << right_size << "  second_len=" << second_len << endl;
    cout << "codeword_len=" << codeword_len << endl;

    assert(codeword_len > t);

    vector<uint64_t> in((size_t)COLUMN_NUMBER * ROW_NUMBER);
    vector<F> compressed((size_t)right_size * ROW_NUMBER, F(0));
    vector<F> codeword((size_t)second_len * ROW_NUMBER, F(0));
    vector<F> rsout(second_len, F(0));

    mt19937_64 rng(12345);
    uniform_int_distribution<uint64_t> dist(0, (uint64_t)1e10);
    for (size_t i = 0; i < (size_t)COLUMN_NUMBER * ROW_NUMBER; i++) {
        in[i] = dist(rng);
    }

    Random_matrix g(COLUMN_NUMBER, right_size, degree, ROW_NUMBER);

    auto t0 = steady_clock::now();

    g.aggregate_batch_left_major(in.data(), compressed.data());
    auto t1 = steady_clock::now();

    cout << "LDPC aggregate time = "
         << duration_cast<milliseconds>(t1 - t0).count() << " ms" << endl;

    for (int b = 0; b < ROW_NUMBER; b++) {
        vector<F> tmp(right_size, F(0));
        for (int i = 0; i < right_size; i++) {
            tmp[i] = compressed[(size_t)i * ROW_NUMBER + b];
        }

        inner_code.encode(tmp.data(), rsout.data());

        for (int i = 0; i < second_len; i++) {
            codeword[(size_t)i * ROW_NUMBER + b] = rsout[i];
        }
    }
    auto t2 = steady_clock::now();

    cout << "Brakedown inner encode time = "
         << duration_cast<milliseconds>(t2 - t1).count() << " ms" << endl;

    Hash* mroot = new Hash[codeword_len];
    for (int i = 0; i < COLUMN_NUMBER; i++) {
        mroot[i] = sha256(
            reinterpret_cast<const uint8_t*>(&in[(size_t)i * ROW_NUMBER]),
            sizeof(uint64_t) * ROW_NUMBER
        );
    }
    for (int i = 0; i < second_len; i++) {
        mroot[i + COLUMN_NUMBER] = sha256(
            reinterpret_cast<const uint8_t*>(&codeword[(size_t)i * ROW_NUMBER]),
            sizeof(F) * ROW_NUMBER
        );
    }

    Hash ROOT = sha256(
        reinterpret_cast<const uint8_t*>(mroot),
        sizeof(Hash) * codeword_len
    );
    (void)ROOT;

    auto t3 = steady_clock::now();

    cout << "Hash time = "
         << duration_cast<milliseconds>(t3 - t2).count() << " ms" << endl;
    cout << "Total Commit time = "
         << duration_cast<milliseconds>(t3 - t0).count() << " ms" << endl;

    F* r  = new F[ROW_NUMBER];
    F* r2 = new F[ROW_NUMBER];
    for (int b = 0; b < ROW_NUMBER; b++) {
        r[b].setByCSPRNG();
        r2[b].setByCSPRNG();
    }

    F* yr = new F[codeword_len];
    F* y1 = new F[codeword_len];
    for (int i = 0; i < codeword_len; i++) {
        yr[i] = 0;
        y1[i] = 0;
    }

    auto t5 = steady_clock::now();
    for (int b = 0; b < ROW_NUMBER; b++) {
        for (int i = 0; i < COLUMN_NUMBER; i++) {
            yr[i] += r[b]  * F(in[(size_t)i * ROW_NUMBER + b]);
            y1[i] += r2[b] * F(in[(size_t)i * ROW_NUMBER + b]);
        }
    }
    auto t6 = steady_clock::now();
    cout << "Open time = "
         << duration_cast<milliseconds>(t6 - t5).count() << " ms" << endl;

    auto t8 = steady_clock::now();

    vector<F> tmp_v(right_size, F(0));
    g.aggregate_Fr(yr, tmp_v.data());

    vector<F> inner_cw(second_len, F(0));
    inner_code.encode(tmp_v.data(), inner_cw.data());
    for (int i = 0; i < second_len; i++) {
        yr[COLUMN_NUMBER + i] = inner_cw[i];
    }

    vector<int> I = sample_positions(codeword_len, t);
    for (int col : I) {
        F targ = yr[col];
        F dot = 0;

        if (col >= COLUMN_NUMBER) {
            for (int j = 0; j < ROW_NUMBER; j++) {
                dot += r[j] * codeword[(size_t)(col - COLUMN_NUMBER) * ROW_NUMBER + j];
            }
        } else {
            for (int j = 0; j < ROW_NUMBER; j++) {
                dot += r[j] * F(in[(size_t)col * ROW_NUMBER + j]);
            }
        }
        assert(dot == targ);
    }

    for (int col = 0; col < codeword_len; col++) {
        Hash now;
        if (col >= COLUMN_NUMBER) {
            now = sha256(
                reinterpret_cast<const uint8_t*>(
                    &codeword[(size_t)(col - COLUMN_NUMBER) * ROW_NUMBER]
                ),
                sizeof(F) * ROW_NUMBER
            );
        } else {
            now = sha256(
                reinterpret_cast<const uint8_t*>(&in[(size_t)col * ROW_NUMBER]),
                sizeof(uint64_t) * ROW_NUMBER
            );
        }
        assert(now == mroot[col]);
    }
    auto t9 = steady_clock::now();

    cout << "Verifier time = "
         << duration_cast<milliseconds>(t9 - t8).count() << " ms" << endl;

    long long PS = 1LL * t * ROW_NUMBER + 2LL * COLUMN_NUMBER + codeword_len;
    cout << "Proof size: " << 32LL * PS / 1048576 << " MB" << endl;

    delete[] mroot;
    delete[] r;
    delete[] r2;
    delete[] yr;
    delete[] y1;
    return 0;
}

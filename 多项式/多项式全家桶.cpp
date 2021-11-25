//#pragma GCC optimize(2)
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

const int N = 3000007;
const int p = 998244353, gg = 3, ig = 332738118, img = 86583718;
const int mod = 998244353;

int qpow(int a, int b) {
    int res = 1;
    while (b) {
        if (b & 1) res = 1ll * res * a % mod;
        a = 1ll * a * a % mod;
        b >>= 1;
    }
    return res;
}

namespace Poly {
#define mul(x, y) (1ll * x * y >= mod ? 1ll * x * y % mod : 1ll * x * y)
#define minus(x, y) (1ll * x - y < 0 ? 1ll * x - y + mod : 1ll * x - y)
#define plus(x, y) (1ll * x + y >= mod ? 1ll * x + y - mod : 1ll * x + y)
#define ck(x) (x >= mod ? x - mod : x)//取模运算太慢了

    typedef vector<int> poly;
    const int G = 3;//根据具体的模数而定，原根可不一定不一样！！！
    //一般模数的原根为 2 3 5 7 10 6
    const int inv_G = qpow(G, mod - 2);
    int RR[N], inv[N];
//    int deer[2][22][N];

    void init(const int t) {//预处理出来NTT里需要的w和wn，砍掉了一个log的时间
//        for (int p = 1; p <= t; ++p) {
//            int buf1 = qpow(G, (mod - 1) / (1 << p));
//            int buf0 = qpow(inv_G, (mod - 1) / (1 << p));
//            deer[0][p][0] = deer[1][p][0] = 1;
//            for (int i = 1; i < (1 << p); ++i) {
//                deer[0][p][i] = 1ll * deer[0][p][i - 1] * buf0 % mod;//逆
//                deer[1][p][i] = 1ll * deer[1][p][i - 1] * buf1 % mod;
//            }
//        }
        inv[1] = 1;
        for (int i = 2; i <= (1 << t); ++i)
            inv[i] = 1ll * inv[mod % i] * (mod - mod / i) % mod;
    }

    int NTT_init(int n) {//快速数论变换预处理
        int limit = 1, L = 0;
        while (limit <= n) limit <<= 1, L++;
        for (int i = 0; i < limit; ++i)
            RR[i] = (RR[i >> 1] >> 1) | ((i & 1) << (L - 1));
        return limit;
    }

//  省空间用
    int deer[2][N];

    void NTT(poly &A, int type, int limit) {//快速数论变换
        A.resize(limit);
        for (int i = 0; i < limit; ++i)
            if (i < RR[i])
                swap(A[i], A[RR[i]]);
        for (int mid = 2, j = 1; mid <= limit; mid <<= 1, ++j) {
            int len = mid >> 1;

//          省空间用
            int buf1 = qpow(G, (mod - 1) / (1 << j));
            int buf0 = qpow(inv_G, (mod - 1) / (1 << j));
            deer[0][0] = deer[1][0] = 1;
            for (int i = 1; i < (1 << j); ++i) {
                deer[0][i] = 1ll * deer[0][i - 1] * buf0 % mod;//逆
                deer[1][i] = 1ll * deer[1][i - 1] * buf1 % mod;
            }

            for (int pos = 0; pos < limit; pos += mid) {
//                int *wn = deer[type][j];
//              省空间用
                int *wn = deer[type];
                for (int i = pos; i < pos + len; ++i, ++wn) {
                    int tmp = 1ll * (*wn) * A[i + len] % mod;
                    A[i + len] = ck(A[i] - tmp + mod);
                    A[i] = ck(A[i] + tmp);
                }
            }
        }
        if (type == 0) {
            for (int i = 0; i < limit; ++i)
                A[i] = 1ll * A[i] * inv[limit] % mod;
        }
    }

    poly poly_mul(poly A, poly B) {//多项式乘法
        int deg = A.size() + B.size() - 1;
        int limit = NTT_init(deg);
        poly C(limit);
        NTT(A, 1, limit);
        NTT(B, 1, limit);
        for (int i = 0; i < limit; ++i)
            C[i] = 1ll * A[i] * B[i] % mod;
        NTT(C, 0, limit);
        C.resize(deg);
        return C;
    }

    poly poly_inv(poly &f, int deg) {//多项式求逆
        if (deg == 1)
            return poly(1, qpow(f[0], mod - 2));

        poly A(f.begin(), f.begin() + deg);
        poly B = poly_inv(f, (deg + 1) >> 1);
        int limit = NTT_init(deg << 1);
        NTT(A, 1, limit), NTT(B, 1, limit);
        for (int i = 0; i < limit; ++i)
            A[i] = B[i] * (2 - 1ll * A[i] * B[i] % mod + mod) % mod;
        NTT(A, 0, limit);
        A.resize(deg);
        return A;
    }

    poly poly_dev(poly f) {//多项式求导
        int n = f.size();
        for (int i = 1; i < n; ++i) f[i - 1] = 1ll * f[i] * i % mod;
        return f.resize(n - 1), f;//f[0] = 0，这里直接扔了,从1开始
    }

    poly poly_idev(poly f) {//多项式求积分
        int n = f.size();
        for (int i = n - 1; i; --i) f[i] = 1ll * f[i - 1] * inv[i] % mod;
        return f[0] = 0, f;
    }

    poly poly_ln(poly f, int deg) {//多项式求对数
        poly A = poly_idev(poly_mul(poly_dev(f), poly_inv(f, deg)));
        return A.resize(deg), A;
    }

    poly poly_exp(poly &f, int deg) {//多项式求指数
        if (deg == 1)
            return poly(1, 1);

        poly B = poly_exp(f, (deg + 1) >> 1);
        B.resize(deg);
        poly lnB = poly_ln(B, deg);
        for (int i = 0; i < deg; ++i)
            lnB[i] = ck(f[i] - lnB[i] + mod);

        int limit = NTT_init(deg << 1);//n -> n^2
        NTT(B, 1, limit), NTT(lnB, 1, limit);
        for (int i = 0; i < limit; ++i)
            B[i] = 1ll * B[i] * (1 + lnB[i]) % mod;
        NTT(B, 0, limit);
        B.resize(deg);
        return B;
    }

    poly poly_sqrt(poly &f, int deg) {//多项式开方
        if (deg == 1) return poly(1, 1);
        poly A(f.begin(), f.begin() + deg);
        poly B = poly_sqrt(f, (deg + 1) >> 1);
        poly IB = poly_inv(B, deg);
        int limit = NTT_init(deg << 1);
        NTT(A, 1, limit), NTT(IB, 1, limit);
        for (int i = 0; i < limit; ++i)
            A[i] = 1ll * A[i] * IB[i] % mod;
        NTT(A, 0, limit);
        for (int i = 0; i < deg; ++i)
            A[i] = 1ll * (A[i] + B[i]) * inv[2] % mod;
        A.resize(deg);
        return A;
    }

    poly poly_pow(poly f, int k) {//多项式快速幂
        if(f.size()==1){
            f[0] = qpow(f[0],k);
            return f;
        }
        f = poly_ln(f, f.size());
        for (auto &x: f) x = 1ll * x * k % mod;
        return poly_exp(f, f.size());
    }

    poly poly_cos(poly f, int deg) {//多项式三角函数（cos）
        poly A(f.begin(), f.begin() + deg);
        poly B(deg), C(deg);
        for (int i = 0; i < deg; ++i)
            A[i] = 1ll * A[i] * img % mod;

        B = poly_exp(A, deg);
        C = poly_inv(B, deg);
        int inv2 = qpow(2, mod - 2);
        for (int i = 0; i < deg; ++i)
            A[i] = 1ll * (1ll * B[i] + C[i]) % mod * inv2 % mod;
        return A;
    }

    poly poly_sin(poly f, int deg) {//多项式三角函数（sin）
        poly A(f.begin(), f.begin() + deg);
        poly B(deg), C(deg);
        for (int i = 0; i < deg; ++i)
            A[i] = 1ll * A[i] * img % mod;

        B = poly_exp(A, deg);
        C = poly_inv(B, deg);
        int inv2i = qpow(img << 1, mod - 2);
        for (int i = 0; i < deg; ++i)
            A[i] = 1ll * (1ll * B[i] - C[i] + mod) % mod * inv2i % mod;
        return A;
    }

    poly poly_arcsin(poly f, int deg) {
        poly A(f.size()), B(f.size()), C(f.size());
        A = poly_dev(f);
        B = poly_mul(f, f);
        for (int i = 0; i < deg; ++i)
            B[i] = minus(mod, B[i]);
        B[0] = plus(B[0], 1);
        C = poly_sqrt(B, deg);
        C = poly_inv(C, deg);
        C = poly_mul(A, C);
        C = poly_idev(C);
        return C;
    }

    poly poly_arctan(poly f, int deg) {
        poly A(f.size()), B(f.size()), C(f.size());
        A = poly_dev(f);
        B = poly_mul(f, f);
        B[0] = plus(B[0], 1);
        C = poly_inv(B, deg);
        C = poly_mul(A, C);
        C = poly_idev(C);
        return C;
    }
}

using namespace Poly;
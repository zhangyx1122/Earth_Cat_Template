#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int mod = 998244353;

void add(int &x, int y) {
    (x += y) >= mod && (x -= mod);
}

void sub(int &x, int y) {
    (x -= y) < 0 && (x += mod);
}

struct FWT {
    int extend(int n) {
        int N = 1;
        for (; N < n; N <<= 1);
        return N;
    }

    void FWTor(std::vector<int> &a, bool rev) {
        int n = a.size();
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1) {
            for (int j = 0; j < n; j += l)
                for (int i = 0; i < m; i++) {
                    if (!rev) add(a[i + j + m], a[i + j]);
                    else sub(a[i + j + m], a[i + j]);
                }
        }
    }

    void FWTand(std::vector<int> &a, bool rev) {
        int n = a.size();
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1) {
            for (int j = 0; j < n; j += l)
                for (int i = 0; i < m; i++) {
                    if (!rev) add(a[i + j], a[i + j + m]);
                    else sub(a[i + j], a[i + j + m]);
                }
        }
    }

    void FWTxor(std::vector<int> &a, bool rev) {
        int n = a.size(), inv2 = (mod + 1) >> 1;
        for (int l = 2, m = 1; l <= n; l <<= 1, m <<= 1) {
            for (int j = 0; j < n; j += l)
                for (int i = 0; i < m; i++) {
                    int x = a[i + j], y = a[i + j + m];
                    if (!rev) {
                        a[i + j] = (x + y) % mod;
                        a[i + j + m] = (x - y + mod) % mod;
                    } else {
                        a[i + j] = 1LL * (x + y) * inv2 % mod;
                        a[i + j + m] = 1LL * (x - y + mod) * inv2 % mod;
                    }
                }
        }
    }

    std::vector<int> Or(std::vector<int> a1, std::vector<int> a2) {
        int n = std::max(a1.size(), a2.size()), N = extend(n);
        a1.resize(N), FWTor(a1, false);
        a2.resize(N), FWTor(a2, false);
        std::vector<int> A(N);
        for (int i = 0; i < N; i++) A[i] = 1LL * a1[i] * a2[i] % mod;
        FWTor(A, true);
        return A;
    }

    std::vector<int> And(std::vector<int> a1, std::vector<int> a2) {
        int n = std::max(a1.size(), a2.size()), N = extend(n);
        a1.resize(N), FWTand(a1, false);
        a2.resize(N), FWTand(a2, false);
        std::vector<int> A(N);
        for (int i = 0; i < N; i++) A[i] = 1LL * a1[i] * a2[i] % mod;
        FWTand(A, true);
        return A;
    }

    std::vector<int> Xor(std::vector<int> a1, std::vector<int> a2) {
        int n = std::max(a1.size(), a2.size()), N = extend(n);
        a1.resize(N), FWTxor(a1, false);
        a2.resize(N), FWTxor(a2, false);
        std::vector<int> A(N);
        for (int i = 0; i < N; i++) A[i] = 1LL * a1[i] * a2[i] % mod;
        FWTxor(A, true);
        return A;
    }
} fwt;

int main() {
    int n;
    scanf("%d", &n);
    n = (1 << n);
    std::vector<int> a1(n), a2(n);
    for (int i = 0; i < n; i++) scanf("%d", &a1[i]);
    for (int i = 0; i < n; i++) scanf("%d", &a2[i]);
    std::vector<int> A;
    A = fwt.Or(a1, a2);
    for (int i = 0; i < n; i++) {
        printf("%d%c", A[i], " \n"[i == n - 1]);
    }
    A = fwt.And(a1, a2);
    for (int i = 0; i < n; i++) {
        printf("%d%c", A[i], " \n"[i == n - 1]);
    }
    A = fwt.Xor(a1, a2);
    for (int i = 0; i < n; i++) {
        printf("%d%c", A[i], " \n"[i == n - 1]);
    }
    return 0;
}
#include <bits/stdc++.h>
using namespace std;
struct SuffixArray {
    static const int MAXN = 1100000;
    char s[MAXN];
    int sa[MAXN], t[MAXN], t1[MAXN], c[MAXN], ra[MAXN], height[MAXN], m;
    inline void init() { memset(this, 0, sizeof(SuffixArray)); }

    inline void get_sa(int n) {
        m = 256;
        int *x = t, *y = t1;
        for (int i = 1; i <= m; i++) c[i] = 0;
        for (int i = 1; i <= n; i++) c[x[i] = s[i]]++;
        for (int i = 1; i <= m; i++) c[i] += c[i - 1];
        for (int i = n; i >= 1; i--) sa[c[x[i]]--] = i;
        for (int k = 1; k <= n; k <<= 1) {
            int p = 0;
            for (int i = n - k + 1; i <= n; i++) y[++p] = i;
            for (int i = 1; i <= n; i++)
                if (sa[i] > k) y[++p] = sa[i] - k;
            for (int i = 1; i <= m; i++) c[i] = 0;
            for (int i = 1; i <= n; i++) c[x[y[i]]]++;
            for (int i = 1; i <= m; i++) c[i] += c[i - 1];
            for (int i = n; i >= 1; i--) sa[c[x[y[i]]]--] = y[i];
            std::swap(x, y);
            p = x[sa[1]] = 1;
            for (int i = 2; i <= n; i++) {
                x[sa[i]] = (y[sa[i - 1]] == y[sa[i]] &&
                            y[sa[i - 1] + k] == y[sa[i] + k])
                               ? p
                               : ++p;
            }
            if (p >= n) break;
            m = p;
        }
    }

    inline void get_height(int n) {
        int i, j, k = 0;
        for (int i = 1; i <= n; i++) ra[sa[i]] = i;
        for (int i = 1; i <= n; i++) {
            if (k) k--;
            int j = sa[ra[i] - 1];
            while (s[i + k] == s[j + k]) k++;
            height[ra[i]] = k;
        }
    }

} SA;   //字符串下标从一开始

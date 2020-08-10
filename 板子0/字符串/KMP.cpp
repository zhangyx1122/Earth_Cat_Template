#include <bits/stdc++.h>
using namespace std;
struct KMP {
    static const int MAXN = 1000010;
    char T[MAXN], P[MAXN];
    int fail[MAXN];
    vector<int> ans;

    void init() { ans.clear(); }

    void get_fail() {
        int m = strlen(P);
        fail[0] = fail[1] = 0;
        for (int i = 1; i < m; i++) {
            int j = fail[i];
            while (j && P[i] != P[j]) j = fail[j];
            fail[i + 1] = (P[i] == P[j] ? j + 1 : 0);
        }
    }

    void find() {
        int n = strlen(T), m = strlen(P);
        get_fail();
        int j = 0;
        for (int i = 0; i < n; i++) {
            while (j && P[j] != T[i]) j = fail[j];
            if (P[j] == T[i]) j++;
            if (j == m) ans.push_back(i - m + 1);
        }
    }
} kmp;  //P为模式串，下标从0开始，输入后直接调用find()

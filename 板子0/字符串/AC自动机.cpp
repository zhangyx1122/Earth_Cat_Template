#include <bits/stdc++.h>
using namespace std;
struct AC {
    static const int maxnode = 200005;
    static const int sigma_size = 26;
    char T[maxnode];
    int ch[maxnode][sigma_size];
    int val[maxnode], fail[maxnode], last[maxnode];
    int sz;
    vector<pair<int, int> > ans;

    void init() {
        sz = 1;
        memset(ch[0], 0, sizeof(ch[0]));
        ans.clear();
    }

    int idx(const char &c) { return c - 'a'; }

    void insert(string s, int v) {
        int u = 0, n = s.length();
        for (int i = 0; i < n; i++) {
            int c = idx(s[i]);
            if (!ch[u][c]) {
                memset(ch[sz], 0, sizeof(ch[sz]));
                val[sz] = 0;
                ch[u][c] = sz++;
            }
            u = ch[u][c];
        }
        val[u] = v;
    }

    void get_fail() {
        queue<int> que;
        fail[0] = 0;
        for (int c = 0; c < sigma_size; c++) {
            int u = ch[0][c];
            if (u) {
                fail[u] = 0;
                que.push(u);
                last[u] = 0;
            }
        }
        while (!que.empty()) {
            int r = que.front();
            que.pop();
            for (int c = 0; c < sigma_size; c++) {
                int u = ch[r][c];
                if (!u) continue;
                que.push(u);
                int v = fail[r];
                while (v && !ch[v][c]) v = fail[v];
                fail[u] = ch[v][c];
                last[u] = val[fail[u]] ? fail[u] : last[fail[u]];
            }
        }
    }

    void print(int j) {
        if (j) {
            ans.push_back(pair<int, int>(j, val[j]));
            print(last[j]);
        }
    }

    void find() {
        int n = strlen(T);
        int j = 0;
        for (int i = 0; i < n; i++) {
            int c = idx(T[i]);
            while (j && !ch[j][c]) j = fail[j];
            j = ch[j][c];
            if (val[j])
                print(j);
            else if (last[j])
                print(last[j]);
        }
    }
} ac;   //字符串下标从0开始

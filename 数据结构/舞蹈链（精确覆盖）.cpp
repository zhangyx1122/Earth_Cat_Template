#include <bits/stdc++.h>
using namespace std;
struct DLX {
    static const int maxn = 1000;     //列的上限
    static const int maxr = 1000;     //解的上限
    static const int maxnode = 5000;  //总结点数上限
    int n, sz;
    int S[maxn];

    int row[maxnode], col[maxnode];
    int L[maxnode], R[maxnode], U[maxnode], D[maxnode];

    int ansd, ans[maxr];

    void init(int n) {
        this->n = n;

        //虚拟节点
        for (int i = 0; i <= n; i++) {
            U[i] = i;
            D[i] = i;
            L[i] = i - 1;
            R[i] = i + 1;
        }
        R[n] = 0;
        L[0] = n;

        sz = n + 1;
        memset(S, 0, sizeof(S));
    }

    void addRow(int r, vector<int> columns) {
        int first = sz;
        for (int i = 0; i < columns.size(); i++) {
            int c = columns[i];
            L[sz] = sz - 1;
            R[sz] = sz + 1;
            D[sz] = c;
            U[sz] = U[c];
            D[U[c]] = sz;
            U[c] = sz;
            row[sz] = r;
            col[sz] = c;
            S[c]++;
            sz++;
        }
        R[sz - 1] = first;
        L[first] = sz - 1;
    }
#define FOR(i, A, s) for (int i = A[s]; i != s; i = A[i])
    void remove(int c) {
        L[R[c]] = L[c];
        R[L[c]] = R[c];
        FOR(i, D, c)
        FOR(j, R, i) {
            U[D[j]] = U[j];
            D[U[j]] = D[j];
            --S[col[j]];
        }
    }

    void restore(int c) {
        FOR(i, U, c)
        FOR(j, L, i) {
            ++S[col[j]];
            U[D[j]] = j;
            D[U[j]] = j;
        }
        L[R[c]] = c;
        R[L[c]] = c;
    }

    // d为递归深度
    bool dfs(int d) {
        if (R[0] == 0) {
            ansd = d;     //找到解
            return true;  //记录解的长度
        }

        //找到S最小的列c
        int c = R[0];
        FOR(i, R, 0) if (S[i] < S[c]) c = i;  //第一个未删除的列

        remove(c);      //删除第c列
        FOR(i, D, c) {  //用结点i所在的行能覆盖的所有其他列
            ans[d] = row[i];
            FOR(j, R, i) remove(col[j]);  //删除结点i所在的能覆的所有其他列
            if (dfs(d + 1)) return true;
            FOR(j, L, i) restore(col[j]);  //恢复结点i所在的行能覆盖的所有其他列
        }
        restore(c);  //恢复第c列

        return false;
    }

    bool solve(vector<int>& v) {
        v.clear();
        if (!dfs(0)) return false;
        for (int i = 0; i < ansd; i++) v.push_back(ans[i]);
        return true;
    }
};
//使用时init初始化，vector中存入r行结点列表用addRow加行，solve(ans)后答案按行的选择在ans中

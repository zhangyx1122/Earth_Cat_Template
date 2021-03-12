struct Trie01 {
    static const int maxnode = 2000005;
    static const int sigma_size = 2;
    int ch[maxnode << 5][sigma_size], val[maxnode << 5];
    int rt[maxnode];
    int sz;

    Trie01() {
        sz = 0;
        memset(ch[0], 0, sizeof(ch[0]));
    }

    void insert(int &now, int pre, int v) {
        now = ++sz;
        for (int i = 30; i >= 0; i--) {
            int k = ((v >> i) & 1);
            ch[now][k] = ++sz;
            ch[now][k ^ 1] = ch[pre][k ^ 1];
            val[ch[now][k]] = val[ch[pre][k]] + 1;
            now = ch[now][k];
            pre = ch[pre][k];
        }
    }
} trie;

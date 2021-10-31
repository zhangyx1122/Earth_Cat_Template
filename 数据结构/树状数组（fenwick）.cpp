template <typename T>
struct fenwick {
    vector<T> fenw;
    int n;

    fenwick(int _n) : n(_n) {
        fenw.resize(n);
    }

    void clear(){
        fenw.clear();
        fenw.resize(n);
    }

    void modify(int x, T v) {
        while (x < n) {
            fenw[x] += v;
            //if(fenw[x]>=mod)fenw[x]-=mod;
            x |= (x + 1);
        }
    }

    T get(int x) {
        T v{};
        while (x >= 0) {
            v += fenw[x];
            //if(v>=mod)v-=mod;
            x = (x & (x + 1)) - 1;
        }
        return v;
    }

    T gets(int l,int r){
        T res=get(r)-get(l-1);
        //if(res<0)res+=mod;
        return res;
    }
};

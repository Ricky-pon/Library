template<typename T, typename F> struct DisjointSparseTable{
    int h, w;
    vector<vector<T>> table;
    T et;
    F f;

    DisjointSparseTable(int n, T et, F f, vector<T> &a): w(n), et(et), f(f){
        h = (n == 1 ? 1 : 32-__builtin_clz(n-1));
        table.resize(h);
        rep(i, h) table[i].resize(w);
        build(a);
    }

    inline void build(vector<T> &a){
        rep(i, h){
            int step = 1<<(i+1);
            for(int l=0, s=(1<<i)-1, t=1<<i, r=(1<<(i+1))-1; l<w; l+=step, s+=step, t+=step, r+=step){
                chmin(s, w-1);
                table[i][s] = a[s];
                rFor(j, s, l) table[i][j] = f(a[j], table[i][j+1]);
                if(s == w-1) break;

                chmin(r, w-1);
                table[i][t] = a[t];
                For(j, t+1, r+1) table[i][j] = f(table[i][j-1], a[j]);
                if(r == w-1) break;
            }
        }
    }

    T query(int l, int r){
        if(l == r) return et;
        --r;
        if(l == r) return table[0][l];
        int t = 31-__builtin_clz(l^r);
        return f(table[t][l], table[t][r]);
    }
};
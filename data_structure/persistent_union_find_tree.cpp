struct PersistentUnionFindTree {
    PersistentArray<int> par;

    PersistentUnionFindTree(int n) {
        vector<int> a(n, -1);
        par.build(a);
    }

    bool is_root(int time, int x) { return par.get(time, x) < 0; }

    int find(int time, int x) {
        if (is_root(time, x)) return x;
        return find(time, par.get(time, x));
    }

    int size(int time, int x) { return -par.get(time, find(time, x)); }

    bool unite(int now, int prv, int x, int y) {
        x = find(prv, x);
        y = find(prv, y);
        if (x == y) {
            par.roots[now] = par.roots[prv];
            return false;
        }
        int sx = size(prv, x), sy = size(prv, y);
        if (sx < sy) swap(x, y), swap(sx, sy);
        par.set(now, prv, x, -sx - sy);
        par.set(now, prv, y, x);
        return true;
    }

    bool same(int time, int x, int y) { return find(time, x) == find(time, y); }
};
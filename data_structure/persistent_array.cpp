template <typename T>
struct PersistentArray {
    struct Node {
        T val;
        Node *ch[2];

        Node() {
            val = {};
            rep(c, 2) ch[c] = nullptr;
        }
    };

    unordered_map<int, Node *> roots;

    PersistentArray() {}

    PersistentArray(vector<T> &a) { build(a); }

    PersistentArray(int n, T val) {
        vector<T> a(n, val);
        build(a);
    }

    void build(vector<T> &a) {
        Node *root = new Node();
        roots[-1] = root;
        rep(i, a.size()) {
            Node *node = root;
            int idx = i;
            while (true) {
                if (idx == 0) {
                    node->val = a[i];
                    break;
                }
                if (!node->ch[idx & 1]) node->ch[idx & 1] = new Node();
                node = node->ch[idx & 1];
                idx >>= 1;
            }
        }
    }

    void set(int now, int prv, int idx, T val) {
        if (roots.find(now) == roots.end()) {
            roots[now] = new Node();
            roots[now]->val = roots[prv]->val;
            rep(c, 2) roots[now]->ch[c] = roots[prv]->ch[c];
        }
        Node *now_node = roots[now], *prv_node = roots[prv];
        while (true) {
            if (idx == 0) {
                now_node->val = val;
                return;
            }
            if (!now_node->ch[idx & 1] ||
                now_node->ch[idx & 1] == prv_node->ch[idx & 1]) {
                now_node->ch[idx & 1] = new Node();
                now_node->ch[idx & 1]->val = prv_node->ch[idx & 1]->val;
                rep(c, 2) now_node->ch[idx & 1]->ch[c] =
                    prv_node->ch[idx & 1]->ch[c];
            }
            now_node = now_node->ch[idx & 1];
            prv_node = prv_node->ch[idx & 1];
            idx >>= 1;
        }
    }

    T get(int time, int idx) {
        Node *node = roots[time];
        while (true) {
            if (idx == 0) return node->val;
            node = node->ch[idx & 1];
            idx >>= 1;
        }
    }
};
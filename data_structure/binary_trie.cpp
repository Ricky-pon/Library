template<typename T, int len=31> struct BinaryTrie{
    struct Node{
        Node *l_ptr, *r_ptr;
        int cnt;
        T lxor;
        Node(): l_ptr(nullptr), r_ptr(nullptr), cnt(0), lxor(0){}
    };
 
    Node *root;
    BinaryTrie(): root(new Node()){}

    int count_child(Node *t){
        return t ? t->cnt : 0;
    }

    void push(Node *t, int dep){
        if(t->lxor >> dep & 1) swap(t->l_ptr, t->r_ptr);
        if(t->l_ptr) t->l_ptr->lxor ^= t->lxor;
        if(t->r_ptr) t->r_ptr->lxor ^= t->lxor;
        t->lxor = 0;
    }

    int size(){
        return count_child(root);
    }
 
    Node *_add(Node *t, T val, int dep){
        if(!t) t = new Node();
        if(dep == -1){
            ++t->cnt;
            return t;
        }
        push(t, dep);
        if(!(val>>dep & 1)) t->l_ptr = _add(t->l_ptr, val, dep-1);
        else t->r_ptr = _add(t->r_ptr, val, dep-1);
        t->cnt = count_child(t->l_ptr) + count_child(t->r_ptr);
        return t;
    }
 
    Node *add(int val){
        return _add(root, val, len-1);
    }
 
    Node *_del(Node *t, int val, int dep){
        if(count_child(t) == 0) return t;
        if(dep == -1){
            --t->cnt;
            return t;
        }
        push(t, dep);
        if(!(val>>dep & 1)) t->l_ptr = _del(t->l_ptr, val, dep-1);
        else t->r_ptr = _del(t->r_ptr, val, dep-1);
        t->cnt = count_child(t->l_ptr) + count_child(t->r_ptr);
        return t;
    }
 
    Node *del(int val){
        return _del(root, val, len-1);
    }

    int count(T x){
        auto t = root;
        rrep(i, len){
            if(!t) return 0;
            push(t, i);
            if(!(x>>i & 1)) t = t->l_ptr;
            else t = t->r_ptr;
        }
        return count_child(t);
    }
 
    T get_kth(int K){
        ++K;
        T ret = 0;
        auto t = root;
        rrep(i, len){
            push(t, i);
            if(count_child(t->l_ptr) >= K){
                t = t->l_ptr;
                ret <<= 1;
            }
            else{
                K -= count_child(t->l_ptr);
                t = t->r_ptr;
                ret <<= 1;
                ++ret;
            }
        }
        return ret;
    }

    void set_xor(T x){
        root->lxor ^= x;
    }
};
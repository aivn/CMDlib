#ifndef TENSOR_HPP
#define TENSOR_HPP

#include <vector>
#include <iterator>
#include <algorithm>

int pos2off(const std::vector<int> & v, const std::vector<int> & mul, int D){
    int res = 0;
    for(int i = 0; i < D; i++){
        res += v[i]*mul[i];
    }
    return res;
}

std::vector<int> off2pos(int off, const std::vector<int> & size, int D){
    std::vector<int> res(D);
    for(int i = 0; i < D; i++){
        int num = off%size[i];
        res[i] = num;
        off-=num;
    }
    return res;
}

//вектора glue1, glue2 должны быть сортированы по возрастанию
std::vector<int> glue_vectors(const std::vector<int> & glue1, const std::vector<int> & v1, const std::vector<int> & glue2, const std::vector<int> & v2){
    int n1 = v1.size(); int n2 = v2.size(); int m = glue1.size();
    int g1 = glue1.size(); int g2 = glue2.size();
    int D = n1+n2+m; std::vector<int> res;
    int count1 = 0; int count2 = 0;
    for(int i = 0; i < n1; i++){
        bool b;
        if(count1 < g1){
            if(i != glue1[count1]){
                b = false;
            }else{
                b = true;
            }
        }else{
            b = false;
        }

        if(!b){
            res.push_back(v1[i]);
        }else{
            count1++;
        }
    }
    for(int i = 0; i < n2; i++){
        bool b;
        if(count2 < g2){
            if(i != glue2[count2]){
                b = false;
            }else{
                b = true;
            }
        }else{
            b = false;
        }

        if(!b){
            res.push_back(v2[i]);
        }else{
            count2++;
        }
    }
    return res;
}

template <typename T>
class Tensor{
    std::vector<T> data;
    int D = 0;
    int sz;
    std::vector<int> size;
    std::vector<int> mul;
    std::vector<int> marks;
public:
    Tensor<T> & operator = (const Tensor<T> & t){
        data = t.data; D = t.D; sz = t.sz; size = t.size; mul = t.mul; marks = t.marks;
        return *this;
    }

    void init(const std::vector<int> & size_){
        D = size_.size();
        size.resize(D); mul.resize(D);
        sz = 1; 
        for(int i = 0; i < D; i++){
            size[i] = size_[i];
            sz *= size[i];
            if(i > 0){
                mul[i] = mul[i-1]*size[i-1];
            }else{
                mul[0] = 1;
            }
        }
        data.resize(sz);
        marks.resize(D);
    }
    T& operator [](const std::vector<int> & v){ return data[pos2off(v, mul, D)]; }
    T operator [](const std::vector<int> & v) const{ return data[pos2off(v, mul, D)]; } 

    T& operator [](int off_){ return data[off_]; }
    T operator [](int off_) const{ return data[off_]; } 

    int get_D() const { return D; }
    std::vector<int> get_size() const { return size; }
    std::vector<int> get_mul() const { return mul; }
    int get_sz() const { 
        return sz; 
    }
    void set_mark(int index, int mark){ marks[index] = mark; }
    std::vector<int> get_marks() const{ return marks; }


    class iterator : public std::iterator<std::input_iterator_tag, T> {
    private:
        int off;
        std::vector<int> v;
        std::vector<int> w;
        bool fixed_indexes = false;
        int D;
        std::vector<int> size;
        std::vector<int> mask;
    public:
        const Tensor<T>* tensor = nullptr;
        const std::vector<int> & get_pos() {
            if(!fixed_indexes){
                return v;
            }else{
                for(int i = 0; i < D; i++){
                    w[mask[i]] = v[i];
                }
                return w; 
            } 
        }

        void init(int off_){ 
            off = off_; 
            size = tensor->get_size(); 
            D = tensor->get_D();
            v = off2pos(off, size, D); 
        }

        void init(const std::vector<int> & fixed_indexes_, int off_, bool is_end){
            fixed_indexes = true;
            off = off_;
            std::vector<int> size_big = tensor->get_size();
            w.resize(tensor->get_D());
            D = 0;
            for(int i = 0, sz_ = fixed_indexes_.size(); i < sz_; i++){
                if(fixed_indexes_[i] == -1){
                    size.push_back(size_big[i]);
                    mask.push_back(i);
                    D++;
                }else{
                    w[i] = fixed_indexes_[i];
                }
            }
            if(is_end){
                off = 1;
                for(int i = 0; i < D; i++){
                    off *= size[i];
                }
            }
            v = off2pos(off, size, D); 
        }
        
        bool operator==(const iterator& other) const { return off == other.off; }
        bool operator!=(const iterator& other) const { return !(*this == other); }
        
        int operator*() const { return off; }
        
        void next_pos(){
            for(int i = 0; i < D; i++){
                if(v[i] < size[i]-1){
                    v[i]++;
                    break;
                }else{
                    v[i] = 0;
                }
            }
            off++;
        }

        iterator& operator++() {
            next_pos();
            return *this;
        }
        
        iterator operator++(int) {
            iterator temp = *this;
            next_pos();
            return temp;
        }

        
    };

    iterator begin() const { 
        iterator it; 
        it.tensor = this; 
        it.init(0); 
        return it;
    }
    iterator end() const { 
        iterator it; 
        it.tensor = this; 
        it.init(get_sz()); 
        return it; 
    }
    iterator begin(const std::vector<int> & fixed_indexes_) const{
        iterator it; 
        it.tensor = this; 
        it.init(fixed_indexes_, 0, false);
        return it; 
    }
    iterator end(const std::vector<int> & fixed_indexes_) const{ 
        iterator it; 
        it.tensor = this; 
        it.init(fixed_indexes_, 0, true); 
        return it; 
    }
};

template <typename T>
inline Tensor<T> tensor_convolution(const Tensor<T> & T1, const std::vector<int> & v1, const Tensor<T> & T2, const std::vector<int> & v2){
    int n = v1.size(); std::vector<int> sz1 = T1.get_size(); std::vector<int> sz2 = T2.get_size();
    int D1 = T1.get_D(); int D2 = T2.get_D(); int D = D1+D2-2*n;
    std::vector<int> v1_sorted = v1; std::vector<int> v2_sorted = v2;
    std::sort(v1_sorted.begin(), v1_sorted.end()); std::sort(v2_sorted.begin(), v2_sorted.end());
    std::vector<int> size = glue_vectors(v1_sorted, sz1, v2_sorted, sz2);
    std::vector<int> new_marks = glue_vectors(v1_sorted, T1.get_marks(), v2_sorted, T2.get_marks());
    
    Tensor<T> res; res.init(size);
    for(int i = 0; i < D; i++){
        res.set_mark(i, new_marks[i]);
    }

    typename Tensor<T>::iterator it1 = T1.begin();
    typename Tensor<T>::iterator it1_end = T1.end();
    while(it1 != it1_end){
        std::vector<int> fixed_indexes(D2);
        for(int i = 0; i < D2; i++){
            fixed_indexes[i] = -1;
        }
        const std::vector<int> & w1 = it1.get_pos();
        for(int i = 0, sz_ = v1.size(); i < sz_; i++){
            fixed_indexes[v2[i]] = w1[v1[i]];
        }
        typename Tensor<T>::iterator it2 = T2.begin(fixed_indexes);
        typename Tensor<T>::iterator it2_end = T2.end(fixed_indexes);
        while(it2 != it2_end){
            const std::vector<int> & w2 = it2.get_pos();
            std::vector<int> w = glue_vectors(v1_sorted, w1, v2_sorted, w2);
            res[w] += T1[w1]*T2[w2];
            it2++;
        }
        it1++;
    }

    return res;
}

template <typename T>
inline Tensor<T> tensor_convolution(const Tensor<T> & t, const std::vector<int> & v1, const std::vector<int> & v2){
    int n{v1.size()}, m{t.get_D()};
    std::vector<int> table1(m-2*n), table2(n), table3(n);
    std::vector<int> size1, size2;
    std::vector<int> t_size = t.get_size();
    for(int i = 0; i < n; i++){
        table2[i] = v1[i];
        table3[i] = v2[i];
        size2.push_back(t_size[v1[i]]);
    }

    std::vector<int> v3;
    for(int i = 0; i < n; i++){
        v3.push_back(v1[i]);
    }
    for(int i = 0; i < n; i++){
        v3.push_back(v2[i]);
    }
    std::vector<int> v3_sorted = v3;
    std::sort(v3_sorted.begin(), v3_sorted.end());
    int count = 0; int count1 = 0;
    for(int i = 0; i < m; i++){
        bool b;
        if(count < 2*n){
            if(i != v3_sorted[count]){
                b = false;
            }else{
                b = true;
            }
        }else{
            b = false;
        }

        if(!b){
            table1[count1++] = i;
            size1.push_back(t_size[v3_sorted[i]]);
        }else{
            count++;
        }
    }
    
    Tensor<T> t1; t1.init(size1); Tensor<T> t2; t2.init(size2);
    std::vector<int> marks = t.get_marks();
    for(int i = 0; i < m-2*n; i++){
        t1.set_mark(i, marks[table1[i]]);
    }

    typename Tensor<T>::iterator it1 = t1.begin(); typename Tensor<T>::iterator it1_end = t1.end();
    while(it1 != it1_end){
        std::vector<int> w1 = it1.get_pos();
        typename Tensor<T>::iterator it2 = t2.begin(); typename Tensor<T>::iterator it2_end = t2.end();
        t1[w1] = 0;
        while(it2 != it2_end){
            std::vector<int> w2 = it2.get_pos();
            std::vector<int> w(m);
            for(int i = 0; i < m-2*n; i++){
                w[table1[i]] = w1[i];
            } 
            for(int i = 0; i < n; i++){
                w[table2[i]] = w2[i]; w[table3[i]] = w2[i];
            }
            t1[w1] += t[w];
            it2++;
        }
        it1++;
    }

    return t1;
}

#endif //TENSOR_HPP
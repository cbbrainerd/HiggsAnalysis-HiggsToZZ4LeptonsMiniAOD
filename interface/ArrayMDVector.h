#include <vector>
#include <iostream>

//Simpler to maintain version of ArrayVector: ArrayMDVector<T,d> creates a d-dimensional vector of T (i.e. ArrayVector<int,3> is basically std::vector<std::vector<std::vector<int> > >)
//cast_to_vector takes in an ArrayMDVector<T,d>* and casts it to the equivalent std::vector<...>* for functions that need to be passed a std::vector pointer

template <class T,int d>
class ArrayMDVector : public std::vector<typename ArrayMDVector<T,d-1>::type> {
    public:
        typedef std::vector<typename ArrayMDVector<T,d-1>::vector_type> vector_type;
        typedef typename ArrayMDVector<T,d-1>::type content_type;
        typedef ArrayMDVector<T,d> type;
        content_type& operator[](std::size_t pos);
};

template<class T>
class ArrayMDVector<T,0> {
    public:
        typedef T vector_type;
        typedef T type;
};

template <class T,int d>
typename ArrayMDVector<T,d>::vector_type* cast_to_vector(ArrayMDVector<T,d>* arrayVector) {
    return (typename ArrayMDVector<T,d>::vector_type*) arrayVector;
}

template <class T, int d>
typename ArrayMDVector<T,d>::content_type& ArrayMDVector<T,d>::operator[](std::size_t pos) {
    std::size_t current_length=this->size();
    if(pos==current_length) {
        this->push_back(typename ArrayMDVector<T,d>::content_type());
    }
    else if(pos>current_length) {
        std::cout << "Warning: assigned past end of ArrayMDVector. Previous length: " << current_length << " Position: " << pos << ". This is probably a bug.\n";
        for(std::size_t i=0;i<=pos-current_length;++i) {
            this->push_back(typename ArrayMDVector<T,d>::content_type());
        }
    }
    return std::vector<typename ArrayMDVector<T,d>::content_type>::operator[](pos);
}

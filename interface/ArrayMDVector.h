#include <vector>
#include <iostream>

/* ArrayMDVector.h
 * Author: Christopher Brainerd - University of California, Davis
 * ArrayMDVector<T,d,dv> creates a d-dimensional vector of T (i.e. ArrayVector<int,3> is basically std::vector<std::vector<std::vector<int> > >)
 * cast_to_vector takes in an ArrayMDVector<T,d,dv>* and casts it to the equivalent std::vector<...>* for functions that need to be passed a std::vector pointer
 * Unlike a std::vector, ArrayMDVector can be assigned to as though it were an C-style array
 * Adds default-constructed elements with a value of default_value if assigned past the end. (i.e. ArrayVector<int,1> foo; foo[3]=1; will give the same result as ArrayVector<int,1> foo {0,0,0,1};) Also warns in this case.
 * Usage: 
 *   Variable declaration: ArrayVectorMD<int,1> foo;
 *   Branch declaration: tree->Branch("foo",cast_to_vector(&foo));
 *   Adding elements: for(int i=0;i<bar;++i) foo[i]=i;
 */

template <class T,int d,int dv=0>
class ArrayMDVector : public std::vector<typename ArrayMDVector<T,d-1,dv>::type> {
    public:
        typedef std::vector<typename ArrayMDVector<T,d-1,dv>::vector_type> vector_type;
        typedef typename ArrayMDVector<T,d-1,dv>::type content_type;
        typedef ArrayMDVector<T,d,dv> type;
        typename ArrayMDVector<T,d,dv>::reference operator[](std::size_t pos);
        content_type default_val() { return default_value((ArrayMDVector<T,d,dv>*)nullptr); }
};

template<class T,int dv>
class ArrayMDVector<T,0,dv> {
    public:
        typedef T vector_type;
        typedef T type;
};

//When operator[] is called out of bounds, construct defaults up to that point
//For a vector type, this is just the default constructor for the vector
//For a floating point type, this is -999 for backwards compatibility
//For a bool, this is false for backwards compatibility

template <class T,int d,int dv>
inline typename ArrayMDVector<T,d,dv>::content_type default_value(ArrayMDVector<T,d,dv>*) { return typename ArrayMDVector<T,d,dv>::content_type(); }

template <class T,int dv>
inline T default_value(ArrayMDVector<T,1,dv>*) { return T(dv); }

template <class T,int d,int dv>
typename ArrayMDVector<T,d,dv>::vector_type* cast_to_vector(ArrayMDVector<T,d,dv>* arrayVector) {
    return (typename ArrayMDVector<T,d,dv>::vector_type*) arrayVector;
}

template <class T, int d,int dv>
typename ArrayMDVector<T,d,dv>::reference ArrayMDVector<T,d,dv>::operator[](std::size_t pos) {
    std::size_t current_length=this->size();
    if(pos==current_length) {
        this->push_back(ArrayMDVector<T,d,dv>::default_val());
    }
    else if(pos>current_length) {
        std::cout << "Warning: assigned past end of ArrayMDVector. Previous length: " << current_length << " Position: " << pos << ". This is probably a bug.\n";
        for(std::size_t i=0;i<=pos-current_length;++i) {
            this->push_back(ArrayMDVector<T,d,dv>::default_val());
        }
    }
    return std::vector<typename ArrayMDVector<T,d,dv>::content_type>::operator[](pos);
}

#include <vector>
#include <iostream>

template <class T>
class ArrayVector : public std::vector<T> {
public:
    T& operator[](std::size_t pos);
};

template <class T>
T& ArrayVector<T>::operator[](std::size_t pos) {
    std::size_t current_length=this->size();
    if(pos==current_length) {
        this->push_back(T());
    }
    if(pos>current_length) {
        for(std::size_t i=0;i<pos-current_length+1;++i) {
            this->push_back(T());
        }
        std::cout << "Warning: vector assigned out of order. Position: " << pos << " Previous length: " << current_length << std::endl;
    }
    return std::vector<T>::operator[](pos);
};

template <class T>
class ArrayVectorDebug : public ArrayVector<T> {
    public:
        T&operator[](std::size_t pos) {
            std::cout << "pos" << pos << "\n";
            return ArrayVector<T>::operator[](pos);
        }
        void clear() {
            std::cout << "Cleared.\n";
            ArrayVector<T>::clear();
        }
}


#include <vector>

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
    return std::vector<T>::operator[](pos);
}

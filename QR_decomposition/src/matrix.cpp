#include "matrix.h"
#include "common.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <utility>

template <typename T>
TMatrix<T>::TMatrix(size_t size1, size_t size2)
    : Size1{size1}
    , Size2{size2}
    {
        Data = new T[Size1 * Size2];
    }

template <typename T>
TMatrix<T>::TMatrix(size_t size)
    : TMatrix(size, size)
    {
    }

template <typename T>
TMatrix<T>::TMatrix(const TMatrix& other) 
    : TMatrix<T>(other.Size1, other.Size2)
    {
        for (size_t i = 0; i < Size1; ++i) {
            for (size_t j = 0; j < Size2; ++j) {
                Data[i * Size2 + j] = other.Data[i * Size2 + j];
            }
        }
    }

template <typename T>
TMatrix<T>::TMatrix(TMatrix&& other) noexcept 
    : Size1{other.GetSize1()}
    , Size2{other.GetSize2()}
    {
        Data = other.Data;
        other.Data = nullptr;
    }

template <typename T>
TMatrix<T> TMatrix<T>::operator =(const TMatrix& other) {
    Size1 = other.Size1;
    Size2 = other.Size2;
    Data = new T[Size1 * Size2];
    for (size_t i = 0; i < Size1; ++i) {
        for (size_t j = 0; j < Size2; ++j) {
            Data[i * Size2 + j] = other.Data[i * Size2 + j];
        }
    }
    return *this;
}

template <typename T>
TMatrix<T> TMatrix<T>::operator =(TMatrix&& other) noexcept {
    Size1 = other.Size1;
    Size2 = other.Size2;
    Data = other.Data;
    other.Data = nullptr;
    return *this;
}

template <typename T>
T* TMatrix<T>::GetData() const noexcept {
    return Data;
}

template <typename T>
size_t TMatrix<T>::GetSize1() const noexcept {
    return Size1;
}

template <typename T>
size_t TMatrix<T>::GetSize2() const noexcept {
    return Size2;
}

template <typename T>
bool TMatrix<T>::operator !() const noexcept {
    return !Data;
}

template <typename T>
T& TMatrix<T>::operator ()(size_t row, size_t column) const {
    if (!Data) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (row > Size1 || column > Size2) {
        throw(std::out_of_range("Cannot get matrix element."));
    }
    return Data[row * Size2 + column];
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator +=(const TMatrix<T>& other) {
    if (!Data || !other.Data) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (other.GetSize1() != Size1 || other.GetSize2() != Size2) {
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply +=."));
    }
    for (size_t i = 0; i < Size1; ++i) {
        for (size_t j = 0; j < Size2; ++j) {
            Data[i * Size1 + j] += other(i, j);
        }
    }
    return *this;
}

template <typename T>
TMatrix<T>& TMatrix<T>::operator -=(const TMatrix<T>& other) {
    if (!Data || !other.Data) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (other.GetSize1() != Size1 || other.GetSize2() != Size2) {
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply -=."));
    }
    for (size_t i = 0; i < Size1; ++i) {
        for (size_t j = 0; j < Size2; ++j) {
            Data[i * Size1 + j] -= other(i, j);
        }
    }
    return *this;
}

template <typename T> template <typename DType>
TMatrix<T>& TMatrix<T>::operator *=(const DType coeff) {
    if (!Data) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    for (size_t i = 0; i < Size1; ++i) {
        for (size_t j = 0; j < Size2; ++j) {
            Data[i * Size1 + j] *= coeff;
        }
    }
    return *this;
}

template <typename T>
void TMatrix<T>::Nullify() {
    for (size_t i = 0; i < Size1; ++i) {
        for (size_t j = 0; j < Size2; ++j) {
            Data[i * Size1 + j] = 0;
        }
    }
}

template <typename T>
TMatrix<T> TMatrix<T>::Transpose() {
    TMatrix<T> res(Size2, Size1);
    for (size_t i = 0; i < Size1; ++i) {
        for (size_t j = 0; j < Size2; ++j) {
            res.Data[j * Size1 + i] = Data[i * Size2 + j];
        }
    }
    return res;
}

template <typename T>
bool TMatrix<T>::IsTriangular(const TMatrix<T>& a, ETriangularType type) {
    if (!a) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    bool isUpper = type == ETriangularType::Upper;
    bool flag = true;
    for (size_t i = 0; i < a.GetSize1(); ++i) {
        for (size_t j = (isUpper ? 0 : i + 1); j < (isUpper ? i : a.GetSize2()); ++j) {
            flag &= RoughEq(a(i, j), 0.0);
        }
    }
    return flag;
}


template <typename T>
double TMatrix<T>::ColumnNorm2(const TMatrix& a, size_t colNum) {
    if (!a) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (colNum > a.GetSize2()) {
        throw(std::out_of_range("Cannot get matrix column."));
    }
    double norm = 0;
    for (size_t i = 0; i < a.GetSize1(); ++i) {
        norm += a(i, colNum) * a(i, colNum);
    }
    return std::sqrt(norm);
}

template <typename T>
TMatrix<T> TMatrix<T>::Prod(const TMatrix<T>& a, const TMatrix<T>& b) {
    if (!a || !b) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (a.GetSize2() != b.GetSize1()) {
        throw(std::invalid_argument("Matrices have no incorrect size. Cannot apply Prod"));
    }
    TMatrix<T> res(a.GetSize1(), b.GetSize2());
    for (size_t i = 0; i < a.GetSize1(); ++i) {
        for (size_t j = 0; j < b.GetSize2(); ++j) {
            res(i, j) = 0;
            for (size_t k = 0; k < a.GetSize2(); ++k) {
                res(i, j) += a(i, k) * b(k, j);
            }
        }
    }
    return res;
}

template <typename T>
TMatrix<T> TMatrix<T>::CreateIdentityMatrix(const size_t n) {
    TMatrix<T> I(n);
    I.Nullify();
    for (size_t i = 0; i < n; ++i) {
        I(i, i) = 1;
    }
    return I;
}

template <typename T>
TMatrix<T> TMatrix<T>::CreateRandom(const size_t size1, const size_t size2) {
    std::srand(std::time(nullptr));
    TMatrix<double> res(size1, size2);
    for (size_t i = 0; i < size1; ++i) {
        for (size_t j = 0; j < size2; ++j) {
            res(i, j) = std::rand() % (size1 * size2);
        }
    }
    return res;
}

template <typename T>
T TMatrix<T>::InnerProd(const TMatrix<T>& a, size_t a_column, const TMatrix<T>& b, size_t b_column) {
    if (!a || !b) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (a.GetSize1() != b.GetSize1()) {
        throw(std::invalid_argument("Inappropriate matrix sizes."));
    }
    T result{};
    for (size_t i = 0; i < a.GetSize1(); ++i) {
        result += a(i, a_column) * b(i, b_column);
    }
    return result;
}

template <typename T>
TMatrix<T>::~TMatrix() {
    if (Data != nullptr) {
        delete[] Data;
    }
}

template <typename T>
std::ostream& operator <<(std::ostream& out, const TMatrix<T>& matrix) {
    if (!matrix) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    for (size_t i = 0; i < matrix.GetSize1(); ++i) {
        for (size_t j = 0; j < matrix.GetSize2(); ++j) {
            out << matrix(i, j) << ' ';
        }
    }
    out <<std::endl;
    return out;
}

template <typename T>
bool operator ==(const TMatrix<T>& a, const TMatrix<T>& b) {
    if (!a || !b) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (a.GetSize1() != b.GetSize2() || a.GetSize2() != b.GetSize2()) {
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply +."));
    }
    for (size_t i = 0; i < a.GetSize1(); ++i) {
        for (size_t j = 0; j < a.GetSize2(); ++j) {
            if (!RoughEq(a(i, j), b(i, j))) {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
bool operator !=(const TMatrix<T>& a, const TMatrix<T>& b) {
    if (!a || !b) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (a.GetSize1() != b.GetSize2() || a.GetSize2() != b.GetSize2()) {
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply +."));
    }
    return !(a == b);
}

template <typename T>
TMatrix<T> operator+(const TMatrix<T>& a, const TMatrix<T>& b) {
    if (!a || !b) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (a.GetSize1() != b.GetSize2() || a.GetSize2() != b.GetSize2()) {
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply +."));
    }
    TMatrix res(a);
    return res += b;
}

template <typename T>
TMatrix<T> operator -(const TMatrix<T>& a, const TMatrix<T>& b) {
    if (!a || !b) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (a.GetSize1() != b.GetSize2() || a.GetSize2() != b.GetSize2()) {
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply +."));
    }
    TMatrix res(a);
    return res -= b;
}

template <typename T1, typename T2>
TMatrix<T1> operator *(const TMatrix<T1>& a, const T2 coeff) {
    if (!a) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    TMatrix res(a);
    return res *= coeff;
}

template class TMatrix<double>;

template std::ostream& operator <<(std::ostream& out, const TMatrix<double>& matrix);

template bool operator ==(const TMatrix<double>& a, const TMatrix<double>& b);

template bool operator !=(const TMatrix<double>& a, const TMatrix<double>& b);

template TMatrix<double> operator +(const TMatrix<double>& a, const TMatrix<double>& b);

template TMatrix<double> operator -(const TMatrix<double>& a, const TMatrix<double>& b);

template TMatrix<double> operator *(const TMatrix<double>& a, const double coeff);
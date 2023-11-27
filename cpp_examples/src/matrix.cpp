#include "common.h"
#include "matrix.h"
#include "vector.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <future>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

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
TMatrix<T>& TMatrix<T>::operator =(const TMatrix& other) {
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
TMatrix<T>& TMatrix<T>::operator =(TMatrix&& other) noexcept {
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
            Data[i * Size2 + j] += other(i, j);
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
TMatrix<T>::operator double() const {
    return Norm2(*this);
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
TMatrix<T> TMatrix<T>::Transpose() const {
    TMatrix<T> res(Size2, Size1);
    for (size_t i = 0; i < Size1; ++i) {
        for (size_t j = 0; j < Size2; ++j) {
            res.Data[j * Size1 + i] = Data[i * Size2 + j];
        }
    }
    return res;
}

template <typename T>
TMatrix<T> TMatrix<T>::CreateFromRange(size_t row1, size_t row2, size_t col1, size_t col2) const { // [row1,row2), [col1, col2)
    if (row1 >= row2 || col1 >= col2) {
        throw(std::invalid_argument("Invalid arguments."));
    }
    TMatrix<T> res(row2 - row1, col2 - col1);
    for (size_t i = row1; i < row2; ++i) {
        for (size_t j = col1; j < col2; ++j) {
            try {
                res(i - row1, j - col1) = Data[i * Size2 + j];
            } catch(...) {
                std::cout << "ERROR: " << i - row1 << ' ' << j - col1 << std::endl;
            }
        }
    }
    return res;
}

template <typename T>
void TMatrix<T>::AssignBlock(TMatrix<T>& matrixBlock, size_t row1, size_t row2, size_t col1, size_t col2) { // [row1,row2), [col1, col2)
    if (row1 > row2 || col1 > col2 || !matrixBlock || Size1 < row2 || Size2 < col2) {
        throw(std::invalid_argument("Invalid arguments."));
    }
    for (size_t i = row1; i < row2; ++i) {
        for (size_t j = col1; j < col2; ++j) {
            Data[i * Size2 + j] = matrixBlock(i - row1, j - col1);
        }
    }
}

template <typename T>
void TMatrix<T>::AssignColumn(size_t columnIdx, const TVector<T>& v) {
    for (size_t i = 0; i < Size1; ++i) {
        Data[i * Size2 + columnIdx] = v(i);
    }
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
double TMatrix<T>::Norm2(const TMatrix& a) {
    double res = 0;
    for (size_t i = 0; i < a.GetSize1(); ++i) {
        for (size_t j = 0; j < a.GetSize2(); ++j) {
            res += a(i, j);
        }
    }
    return std::sqrt(res);
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
bool TMatrix<T>::AbleToMultiply(const TMatrix<T>& a, const TMatrix<T>& b) try {
    if (!a || !b) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    if (a.GetSize2() != b.GetSize1()) {
        throw(std::invalid_argument("Matrices have no incorrect size. Cannot apply Prod"));
    }
    return true;
} catch(std::exception& e) {
    std::cout << e.what() << std::endl;
    return false;
}

template <typename T>
TMatrix<TMatrix<T>> TMatrix<T>::MakeBlockMatrix(const TMatrix& a, const std::pair<size_t, size_t>& blockSize) {
    size_t n = a.GetSize1();
    size_t m = a.GetSize2();
    size_t blockSize1 = blockSize.first;
    size_t blockSize2 = blockSize.second;
    size_t blockMatrixSize1 = (n +  blockSize1 - 1) / blockSize1;
    size_t blockMatrixSize2 = (m +  blockSize2 - 1) / blockSize2;
    TMatrix<TMatrix<T>> blockMatrix(blockMatrixSize1, blockMatrixSize2);
    for (size_t j = 0; j < m; j += blockSize2) {
        for (size_t i = 0; i < n; i += blockSize1) {
            blockMatrix(i / blockSize1, j / blockSize2) = a.CreateFromRange(i, std::min(i + blockSize1, n), j, std::min(j + blockSize2, m));
        }
    }
    return blockMatrix;
}

template <typename T>
TMatrix<T> TMatrix<T>::Prod(const TMatrix<T>& a, const TMatrix<T>& b) {
    if (!AbleToMultiply(a, b)) {
        return TMatrix<T>();
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
TMatrix<T> TMatrix<T>::ParallelProd(const TMatrix& a, const TMatrix& b) {
    if (!AbleToMultiply(a, b)) {
        return TMatrix<T>();
    }
    size_t n = a.GetSize1();
    size_t m = a.GetSize2();
    size_t l = b.GetSize2();
    size_t blockSize1 = std::max<size_t>(std::sqrt(n), 1 << 8);
    size_t blockSize2 = std::max<size_t>(std::sqrt(m), 1 << 8);
    size_t blockSize3 = std::max<size_t>(std::sqrt(l), 1 << 8);
    size_t aBlockMatrixSize1 = (n +  blockSize1 - 1) / blockSize1;
    size_t aBlockMatrixSize2;
    size_t bBlockMatrixSize1 = aBlockMatrixSize2 = (m +  blockSize2 - 1) / blockSize2;
    size_t bBlockMatrixSize2 = (l +  blockSize3 - 1) / blockSize3;

    auto aBlockMatrix = MakeBlockMatrix(a, {blockSize1, blockSize2});
    auto bBlockMatrix = MakeBlockMatrix(a, {blockSize2, blockSize3});

    TMatrix<TMatrix<T>> cBlockMatrix(aBlockMatrixSize1, bBlockMatrixSize2);
    for (size_t i = 0; i < aBlockMatrixSize1; ++i) {
        for (size_t j = 0; j < bBlockMatrixSize2; ++j) {
            std::vector<std::future<TMatrix<T>>> multResults;
            for (size_t k = 0; k < aBlockMatrixSize2; ++k) {
                // std::future<int> f2 = std::async(std::launch::async, []{ return 8; });
                if (!cBlockMatrix(i, j)) {
                    cBlockMatrix(i, j) = TMatrix<T>::Prod(aBlockMatrix(i, k), bBlockMatrix(k, j));
                } else {
                    multResults.emplace_back(std::async(std::launch::async, TMatrix<T>::Prod, std::ref(aBlockMatrix(i, k)), std::ref(bBlockMatrix(k, j))));
                    // cBlockMatrix(i, j) += TMatrix<T>::Prod(aBlockMatrix(i, k), bBlockMatrix(k, j));
                }
            }
            for (auto&& multResult : multResults) {
                cBlockMatrix(i, j) += multResult.get();
            }
        }
    }
    return TMatrix<T>::FromBlockMatrix(cBlockMatrix);
}


template <typename T>
TMatrix<T> TMatrix<T>::BlockProd(const TMatrix<T>& a, const TMatrix<T>& b, const std::tuple<size_t, size_t, size_t>& blockSizes) {
    if (!AbleToMultiply(a, b)) {
        return TMatrix<T>();
    }
    size_t n = a.GetSize1();
    size_t m = a.GetSize2();
    size_t l = b.GetSize2();
    size_t blockSize1 = std::get<0>(blockSizes);
    size_t blockSize2 = std::get<1>(blockSizes);
    size_t blockSize3 = std::get<2>(blockSizes);
    size_t aBlockMatrixSize1 = (n +  blockSize1 - 1) / blockSize1;
    size_t aBlockMatrixSize2;
    size_t bBlockMatrixSize1 = aBlockMatrixSize2 = (m +  blockSize2 - 1) / blockSize2;
    size_t bBlockMatrixSize2 = (l +  blockSize3 - 1) / blockSize3;
    auto aBlockMatrix = MakeBlockMatrix(a, {blockSize1, blockSize2});
    auto bBlockMatrix = MakeBlockMatrix(a, {blockSize2, blockSize3});
    
    TMatrix<TMatrix<T>> cBlockMatrix(aBlockMatrixSize1, bBlockMatrixSize2);
    for (size_t i = 0; i < aBlockMatrixSize1; ++i) {
        for (size_t j = 0; j < bBlockMatrixSize2; ++j) {
            for (size_t k = 0; k < aBlockMatrixSize2; ++k) {
                if (!cBlockMatrix(i, j)) {
                    cBlockMatrix(i, j) = TMatrix<T>::Prod(aBlockMatrix(i, k), bBlockMatrix(k, j));
                } else {
                    cBlockMatrix(i, j) += TMatrix<T>::Prod(aBlockMatrix(i, k), bBlockMatrix(k, j));
                }
            }
        }
    }
    return TMatrix<T>::FromBlockMatrix(cBlockMatrix);
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
    TMatrix<T> res(size1, size2);
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
TMatrix<T> TMatrix<T>::FromBlockMatrix(const TMatrix<TMatrix<T>>& a) {
    if (!a) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    size_t n = 0, m = 0;
    for (size_t i = 0; i < a.GetSize1(); ++i) {
        n += a(i, 0).GetSize1();
    }
    for (size_t j = 0; j < a.GetSize2(); ++j) {
        m += a(0, j).GetSize2();
    }
    TMatrix<T> res(n, m);
    size_t i = 0;
    size_t j = 0;
    for (size_t ii = 0; ii < a.GetSize1(); ++ii) {
        for (size_t jj = 0; jj < a.GetSize2(); ++jj) {
            for (size_t cur_i = i; cur_i - i < a(ii, jj).GetSize1(); ++cur_i) {
                for (size_t cur_j = j; cur_j - j < a(ii, jj).GetSize2(); ++cur_j) {
                    res(cur_i, cur_j) = a(ii, jj)(cur_i - i, cur_j - j);
                }
            }
            j += a(ii, jj).GetSize2();
            if (jj == a.GetSize2() - 1) {
                i += a(ii, jj).GetSize1();
                j = 0;
            }
        }   
    }
    return res;
}


template <typename T> 
TMatrix<T>::~TMatrix() {
    if (Data != nullptr) {
        delete[] Data;
        Data = nullptr;
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
        out << '\n';
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
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply =="));
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
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply !="));
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
        throw(std::invalid_argument("Matrices sizes are different. Cannot apply -."));
    }
    TMatrix res(a);
    return res -= b;
}

template <typename T>
TMatrix<T> operator *(const TMatrix<T>& a, const TMatrix<T>& b) {
    return TMatrix<T>::Prod(a, b); 
}

template <typename T1, typename T2>
TMatrix<T1> operator *(const TMatrix<T1>& a, const T2 coeff) {
    if (!a) {
        throw(std::invalid_argument("Matrix has no Data."));
    }
    TMatrix res(a);
    return res *= coeff;
}

template std::ostream& operator <<(std::ostream& out, const TMatrix<double>& matrix);
template bool operator ==(const TMatrix<double>& a, const TMatrix<double>& b);
template bool operator !=(const TMatrix<double>& a, const TMatrix<double>& b);
template TMatrix<double> operator +(const TMatrix<double>& a, const TMatrix<double>& b);
template TMatrix<double> operator -(const TMatrix<double>& a, const TMatrix<double>& b);
template TMatrix<double> operator *(const TMatrix<double>& a, const TMatrix<double>& b);
template TMatrix<double> operator *(const TMatrix<double>& a, const double coeff);

template std::ostream& operator <<(std::ostream& out, const TMatrix<TMatrix<double>>& matrix);
template bool operator ==(const TMatrix<TMatrix<double>>& a, const TMatrix<TMatrix<double>>& b);
template bool operator !=(const TMatrix<TMatrix<double>>& a, const TMatrix<TMatrix<double>>& b);
template TMatrix<TMatrix<double>> operator +(const TMatrix<TMatrix<double>>& a, const TMatrix<TMatrix<double>>& b);
template TMatrix<TMatrix<double>> operator -(const TMatrix<TMatrix<double>>& a, const TMatrix<TMatrix<double>>& b);
template TMatrix<TMatrix<double>> operator *(const TMatrix<TMatrix<double>>& a, const TMatrix<TMatrix<double>>& b);
template TMatrix<TMatrix<double>> operator *(const TMatrix<TMatrix<double>>& a, const double coeff);

template class TMatrix<double>;
template class TMatrix<TMatrix<double>>;
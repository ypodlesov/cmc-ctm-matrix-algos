#pragma once
#include <iostream>

enum class ETriangularType {
    Upper,
    Lower
};

template <typename T>
class TMatrix {
public:
    TMatrix() = default;
    TMatrix(size_t size1, size_t size2);
    TMatrix(size_t size);
    TMatrix(const TMatrix& other);
    TMatrix(TMatrix& other);
    TMatrix(TMatrix&& other) noexcept;
    TMatrix operator =(const TMatrix& other);
    TMatrix operator =(TMatrix& other);
    TMatrix operator =(TMatrix&& other) noexcept;

    T* GetData() const noexcept;
    size_t GetSize1() const noexcept;
    size_t GetSize2() const noexcept;

    bool operator !() const noexcept;
    T& operator ()(size_t row, size_t column) const;
    TMatrix& operator +=(const TMatrix& other);
    TMatrix& operator -=(const TMatrix& other);
    template <typename DType>
    TMatrix& operator *=(const DType coeff);

    void Nullify();
    TMatrix Transpose();

    static bool IsTriangular(const TMatrix& a, ETriangularType type);
    static double ColumnNorm2(const TMatrix& a, size_t colNum);
    static TMatrix Prod(const TMatrix& a, const TMatrix& b);
    static TMatrix CreateIdentityMatrix(const size_t n);
    static TMatrix CreateRandom(const size_t size1, const size_t size2);


    ~TMatrix();
    
private:
    size_t Size1{};
    size_t Size2{};
    T* Data;
};

template <typename T>
std::ostream& operator <<(std::ostream& out, const TMatrix<T>& matrix);

template <typename T>
bool operator ==(const TMatrix<T>& a, const TMatrix<T>& b);

template <typename T>
bool operator !=(const TMatrix<T>& a, const TMatrix<T>& b);

template <typename T>
TMatrix<T> operator +(const TMatrix<T>& a, const TMatrix<T>& b);

template <typename T>
TMatrix<T> operator -(const TMatrix<T>& a, const TMatrix<T>& b);

template <typename T1, typename T2>
TMatrix<T1> operator *(const TMatrix<T1>& a, const T2 coeff);
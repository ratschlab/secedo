#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <memory>

/**
 * Stores a templatized 2-dimensional matrix, in row-major form.
 * In order to avoid unintentional copying, the class is not assignable or copy-constructable, but
 * it is copy-moveable and move-constructible. You can always make an explicit copy of the class by
 * calling the #copy() method.
 * @tparam T the type of the matrix elements
 */
template <typename T>
class Mat {
  private:
    /**
     * Pointer to a continuous space of memory storing the elements for the Matrix.
     */
    T *el;
    /**
     * If the matrix allocated its own elements, this points to the allocated elements, if the
     * matrix was constructed with an existing memory block, this is nullptr.
     */
    std::unique_ptr<T[]> ownEl;

    /**
     * The size of the matrix (rows, cols).
     */
    uint32_t r, c;


    uint32_t min_idx() const {
        uint32_t result = 0;
        for (uint32_t i = 1; i < r * c; ++i) {
            if (el[i] < el[result]) {
                result = i;
            }
        }
        return result;
    }

    uint32_t max_idx() const {
        uint32_t result = 0;
        for (uint32_t i = 1; i < r * c; ++i) {
            if (el[i] > el[result]) {
                result = i;
            }
        }
        return result;
    }

  public:
    /**
     * Default constructor, sets the matrix to its default state (0 size and null elements).
     */
    Mat() : el(nullptr), ownEl(nullptr), r(0), c(0) {}

    /*
     * Creates a matrix with the given size and elements. The elements are not copied and ownership
     * of the elements remains with the caller.
     */
    Mat(const uint32_t r, const uint32_t c, T *const el) : el(el), ownEl(nullptr), r(r), c(c) {}

    /**
     * Constructs a Mat object from an initializer list.
     * @param r number of rows
     * @param c number of cols
     * @param initEl the elements to copy into the matrix, its size must be r*c
     */
    Mat(const uint32_t r, const uint32_t c, const std::initializer_list<T> &initEl)
        : ownEl(nullptr), r(r), c(c) {
        assert(r * c == initEl.size());
        ownEl = std::unique_ptr<T[]>(new T[r * c]);
        el = ownEl.get();
        std::copy(initEl.begin(), initEl.end(), el);
    }

    /**
     * Creates a matrix with the given size. The matrix will own and manage its elements.
     * @param r number of rows
     * @param c number of cols
     */
    Mat(const uint32_t r, const uint32_t c) : r(r), c(c) {
        ownEl = std::unique_ptr<T[]>(new T[r * c]);
        el = ownEl.get();
    }

    uint32_t rows() const { return r; }
    uint32_t cols() const { return c; }

    /**
     * Delete the copy-constructor. Copy-ing needs to be done explicitly via #copy()
     */
    Mat(const Mat &rhs) = delete;

    /**
     * Move constructor. The elements of the parameter are pilfered and the parameter matrix is set
     * to its default state.
     */
    Mat(Mat &&rhs) : r(rhs.r), c(rhs.c) {
        ownEl = std::move(rhs.ownEl);
        el = rhs.el;
        rhs.r = rhs.c = 0;
        rhs.el = nullptr;
    }

    /**
     * Returns true if the matrix is empty.
     */
    bool empty() const { return r == 0 && c == 0; }

    /**
     * Returns the given element of the matrix.
     */
    T &operator()(uint32_t row, uint32_t col) {
        assert(row < r && col < c && "Invalid matrix index");
        return el[row * c + col];
    }

    /**
     * Returns the given element of the matrix.
     */
    const T &operator()(uint32_t row, uint32_t col) const {
        assert(row < r && col < c && "Invalid matrix index");
        return el[row * c + col];
    }

    /**
     * Returns the sum of two matrices of the same shape
     */
    Mat<T> operator+(Mat<T> other) const {
        assert(r == other.r && c == other.c && "Incompatible matrix shapes for addition");
        Mat<T> result(r, c);
        for (uint32_t i = 0; i < r * c; ++i) {
            result.el[i] = el[i] + other.el[i];
        }
        return result;
    }

    /**
     * Increments the current matrix with elements from the given matrix
     */
    void operator+=(const Mat<T> &other) const {
        assert(r == other.r && c == other.c && "Incompatible matrix shapes for addition");
        for (uint32_t i = 0; i < r * c; ++i) {
            el[i] += other.el[i];
        }
    }

    /**
     * Returns the difference of two matrices of the same shape
     */
    Mat<T> operator-(Mat<T> other) const {
        assert(r == other.r && c == other.c && "Incompatible matrix shapes for subtraction");
        Mat<T> result(r, c);
        for (uint32_t i = 0; i < r * c; ++i) {
            result.el[i] = el[i] - other.el[i];
        }
        return result;
    }

    /** Subtracts #other from the current matrix */
    void operator-=(const Mat<T> &other) {
        assert(r == other.r && c == other.c && "Incompatible matrix shapes for subtraction");
        Mat<T> result(r, c);
        for (uint32_t i = 0; i < r * c; ++i) {
            el[i] -= other.el[i];
        }
    }

    /** Multiply *this with a scalar and return the newly obtained matrix. */
    Mat<T> operator*(const T scalar) const {
        Mat<T> result(r, c);
        for (uint32_t i = 0; i < r * c; ++i) {
            result.el[i] = el[i] * scalar;
        }
        return result;
    }

    /** Multiply with *scalar in-place */
    Mat<T> &operator*=(const T scalar) {
        for (uint32_t i = 0; i < r * c; ++i) {
            el[i] *= scalar;
        }
        return *this;
    }

    /** Add a scalar to *this and return the newly obtained matrix. */
    Mat<T> operator+(const T scalar) const {
        Mat<T> result(r, c);
        for (uint32_t i = 0; i < r * c; ++i) {
            result.el[i] = el[i] + scalar;
        }
        return result;
    }

    /** Add scalar in-place */
    Mat<T> &operator+=(const T scalar) {
        for (uint32_t i = 0; i < r * c; ++i) {
            el[i] += scalar;
        }
        return *this;
    }

    /**
     * Returns the element at the specified index, as stored internally. Eg mat[index] =
     * mat(index/tCol, index%tCol)
     */
    const T &operator[](uint32_t index) const { return el[index]; }
    /**
     * Returns the modifiable element at the specified index, as stored internally. Eg mat[index] =
     * mat(index/tCol, index%tCol)
     */
    T &operator[](uint32_t index) { return el[index]; }

    /**
     * Returns a pointer to the given row of the matrix. Check for boundaries in debug mode.
     */
    T *row(uint32_t row) {
        assert(row < r && "Invalid matrix row");
        return &el[row * c];
    }

    /**
     * Returns a pointer to the given row of the matrix. Check for boundaries in debug mode.
     */
    const T *row(uint32_t row) const {
        assert(row < r && "Invalid matrix row");
        return &el[row * c];
    }

    /**
     * Returns a pointer to the matrix' data
     */
    T *data() { return el; }
    /**
     * Returns a pointer to the matrix' data
     */
    const T *data() const { return el; }

    static uint32_t elemSize() { return sizeof(T); }

    /**
     * The default assign operator is deleted.
     */
    Mat &operator=(const Mat &rhs) = delete;

    /**
     * Move-assign will pilfer the argument, reset it to its default state and copy its contents
     * into the current matrix.
     */
    Mat &operator=(Mat &&rhs) {
        if (this != &rhs) {
            r = rhs.r;
            c = rhs.c;
            ownEl = std::move(rhs.ownEl);
            el = rhs.el;
            rhs.r = 0;
            rhs.c = 0;
            rhs.el = nullptr;
        }
        return *this;
    }

    /**
     * Returns a (deep) copy of the current matrix.
     */
    Mat copy() const {
        Mat result(r, c);
        std::copy(el, el + r * c, result.el);
        return result;
    }

    /**
     * Copies the given row of the current matrix, to the given row of the given matrix
     */
    void copyRow(uint32_t rrc, uint32_t rowDest, Mat<T> *const dest) const {
        std::copy(row(rrc), row(rrc) + c, dest->row(rowDest));
    }

    /**
     * Gets a rectangular block from a matrix of size tRow1xtCol1, starting with point
     * (startRow,startCol). The new matrix is a deep copy of the block in the old matrix.
     *
     * @param startRow first row of the resulting block, such that 0 <= startRow and startRow + nRow
     * <= rows()
     * @param startCol first column of the resulting block, such that 0 <=startCol and startCol +
     * nCol <= cols()
     * @param nr number of rows in the result
     * @param nc number of columns in the result
     */
    Mat block(uint32_t startRow, uint32_t startCol, uint32_t nr, uint32_t nc) const {
        assert(startRow + nr <= r && startCol + nc <= c && "block to extract is out of range");
        Mat result(nr, nc);
        T *ptrDest = result.data();
        const T *ptrSrc = &(*this)(startRow, startCol);
        for (uint32_t row = 0; row < nr; ++row) {
            std::copy(ptrSrc, ptrSrc + nc, ptrDest);
            ptrSrc += c;
            ptrDest += nc;
        }
        return result;
    }

    /**
     * Gets a range of r from the given matrix, starting with startRow and ending with endRow
     * (exclusive).
     * @param startRow first row of the resulting block, such that 0 <= startRow  <= r
     * @param endRow last row of the resulting block, such that 0 <=startRow <= endRow <= r
     */
    Mat rowRange(uint32_t startRow, uint32_t endRow) const {
        assert(startRow <= endRow && endRow <= r && "r to extract are out of range");
        Mat result(endRow - startRow, c);
        T *ptrDest = result.data();
        const T *ptrSrc = &(*this)(startRow, 0);
        std::copy(ptrSrc, ptrSrc + result.r * result.c, ptrDest);
        return result;
    }

    Mat colRange(uint32_t startCol, uint32_t endCol) const {
        assert(startCol <= endCol && endCol <= c && "r to extract are out of range");
        Mat result(r, endCol - startCol);
        T *ptrDest = result.data();
        const T *ptrSrc = &(this->data()[startCol]);
        for (uint32_t row = 0; row < r; ++row) {
            std::copy(ptrSrc, ptrSrc + result.c, ptrDest);
            ptrDest += result.c;
            ptrSrc += c;
        }
        return result;
    }

    /**
     * Returns the transpose of the current matrix.
     */
    Mat transpose() const {
        Mat result(c, r);
        uint32_t idx = 0;
        for (uint32_t row = 0; row < r; ++row) {
            for (uint32_t col = 0; col < c; ++col) {
                result(col, row) = el[idx++];
            }
        }
        return result;
    };

    /*
     * Returns the squared Frobenius norm of the matrix (sum of elements squared), computed as a
     * double.
     */
    double norm2() const {
        double result = 0;
        for (uint32_t i = 0; i < r * c; ++i) {
            result += el[i] * el[i];
        }
        return result;
    }

    /*
     * Returns the Frobenius norm of the matrix (square root of sum of elements squared)
     */
    double norm() const { return std::sqrt(norm2()); }

    /**
     * Returns the position of the largest element in the matrix. If the largest element appears
     * multiple times, the first occurrence is returned.
     */
    std::pair<uint32_t, uint32_t> argMax() {
        uint32_t result = max_idx();
        return { result / c, result % c };
    }

    /**
     * Returns the position of the smallest element in the matrix. If the largest element appears
     * multiple times, the first occurrence is returned.
     */
    std::pair<uint32_t, uint32_t> argMin() const {
        uint32_t result = min_idx();
        return { result / c, result % c };
    }

    /**
     * Inverts each element of the matrix.
     */
    void inv() {
        static_assert(std::is_floating_point_v<T>);
        for (uint32_t i = 0; i < r * c; ++i) {
            el[i] = 1. / el[i];
        }
    }

    /**
     * Exponentiates each element of the matrix.
     */
    void exp() {
        static_assert(std::is_floating_point_v<T>);
        for (uint32_t i = 0; i < r * c; ++i) {
            el[i] = std::exp(el[i]);
        }
    }

    /**
     * Returns the smallest element of the matrix.
     */
    T min() const { return el[min_idx()]; }

    /**
     * Returns the largest element of the matrix.
     */
    T max() const { return el[max_idx()]; }

    /**
     * Set the given row in the matrix to the given value.
     */
    void setRow(uint32_t row, T value) { std::fill_n(el + row * c, c, value); }

    /**
     * Set the given column in the matrix to the given value.
     */
    void setCol(uint32_t col, T value) {
        for (uint32_t row = 0; row < r; ++row) {
            (*this)(row, col) = value;
        }
    }

    /**
     * Checks for matrix element by element equality.
     */
    bool operator==(const Mat<T> &other) const {
        if (this->c != other.c || this->r != other.r) {
            return false;
        }
        if (this == &other || this->el == other.el) {
            return true;
        }
        for (uint32_t i = 0; i < r * c; i++)
            if (el[i] != other.el[i]) {
                return false;
            }
        return true;
    }

    friend std::ostream &operator<<(std::ostream &os, const Mat<T> &mat) {
        if (mat.empty()) {
            os << "[]" << std::endl;
            return os;
        }
        // this operator is only used for debugging so casting is OK
        // this comment is here just so we know we're not stupid :)
        os << "[";
        for (uint32_t i = 0; i < mat.r * mat.c - 1; ++i) {
            os << std::fixed << std::setprecision(5) << (double)mat.el[i];
            if ((i + 1) % mat.c == 0) {
                os << ";\n";
            } else {
                os << ", ";
            }
        }
        os << std::fixed << std::setprecision(5) << (double)mat.el[mat.r * mat.c - 1];
        os << "]\n";
        return os;
    }

    /**
     * Return a matrix with all elements set to value.
     */
    static Mat<T> fill(uint32_t r, uint32_t c, const T &value) {
        Mat<T> result(r, c);
        std::fill_n(result.el, r * c, value);
        return result;
    };

    /**
     * Return a matrix with all elements set to value.
     */
    void fill_diagonal(const T &value) {
        for (uint32_t i = 0; i < std::min(r, c); ++i) {
            (*this)(i, i) = value;
        }
    }

    /**
     * Return a matrix with all elements set to zero.
     */
    static Mat<T> zeros(uint32_t r, uint32_t c) { return Mat::fill(r, c, 0); };

    /**
     * Returns the identity matrix. If the matrix is not square, it initializes the elements on the
     * main diagonal with 1, the rest are zero
     * @return the created identity matrix
     */
    static Mat<T> identity(uint32_t r, uint32_t c) {
        Mat<T> result = Mat<T>::zeros(r, c);
        for (uint32_t i = 0; i < std::min(r, c); ++i) {
            result(i, i) = 1;
        }
        return result;
    };
};

using Matd = Mat<double>;
using Matf = Mat<float>;
using Matu8 = Mat<uint8_t>;
using Mat16u = Mat<uint16_t>;
using Mat32u = Mat<uint32_t>;
using Mat64u = Mat<uint64_t>;
using Matb = Mat<bool>;

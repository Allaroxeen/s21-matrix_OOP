#include "s21_matrix_oop.h"

// Private methods

void S21Matrix::CopyMatrix(const S21Matrix& other) {
  int r = rows_, c = cols_;
  if (other.rows_ < rows_) {
    r = other.rows_;
  }
  if (other.cols_ < cols_) {
    c = other.cols_;
  }
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

void S21Matrix::deallocation() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
}

inline void S21Matrix::InitializeMatrix() {
  matrix_ = new double*[rows_]();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::swaprows(int a, int b) {
  double* temp = matrix_[a];
  matrix_[a] = matrix_[b];
  matrix_[b] = temp;
}

double S21Matrix::MinorMatrixDet(const int a, const int b) {
  S21Matrix minor(rows_ - 1, rows_ - 1);
  for (int i = 0, i2 = 0; ((i < rows_) && (i2 < rows_ - 1)); i++, i2++) {
    if (i == a) {
      i2--;
    } else {
      for (int j = 0, j2 = 0; (j < rows_ && j2 < rows_ - 1); j++, j2++) {
        if (j == b)
          j2--;
        else {
          minor.matrix_[i2][j2] = matrix_[i][j];
        }
      }
    }
  }
  return minor.Determinant();
}

// Constructors/Destructors

S21Matrix::S21Matrix() {
  rows_ = 1;
  cols_ = 1;
  InitializeMatrix();
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows < 1 || cols < 1)
    throw std::runtime_error("Invalid matrix dimension");
  rows_ = rows;
  cols_ = cols;
  InitializeMatrix();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  InitializeMatrix();
  CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix::~S21Matrix() {
  this->deallocation();
  rows_ = 0;
  cols_ = 0;
}

// Setters/Getters

void S21Matrix::SetRows(int rows) {
  if (rows < 1) {
    throw std::runtime_error("Invalid size of rows");
  }
  if (rows_ != rows) {
    S21Matrix temp(rows, cols_);
    temp.CopyMatrix(*this);
    *this = temp;
  }
};
void S21Matrix::SetCols(int cols) {
  if (cols < 1) {
    throw std::runtime_error("Invalid size of columns");
  }
  if (cols_ != cols) {
    S21Matrix temp(rows_, cols);
    temp.CopyMatrix(*this);
    *this = temp;
  };
};
int S21Matrix::GetRows() { return rows_; };
int S21Matrix::GetCols() { return cols_; };

void S21Matrix::PrintMatrix() {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      std::cout << matrix_[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

// Public Methods

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool result = true;
  if (other.rows_ != rows_ || other.cols_ != cols_)
    result = false;
  else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (abs(matrix_[i][j] - other.matrix_[i][j]) > Error_Rate)
          result = false;
      }
    }
  }
  return result;
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::runtime_error("Invalid matrix dimensions for substitute");
  else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] -= other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::runtime_error("Invalid matrix dimensions for summary");
  else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] += other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::runtime_error("Invalid matrix dimensions for multiplication");
  }
  S21Matrix buf(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        buf.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = buf;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix newmatrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      newmatrix.matrix_[j][i] = matrix_[i][j];
    }
  }
  return newmatrix;
}

double S21Matrix::Determinant() {
  double det = 0.0;
  if (rows_ != cols_) {
    throw std::runtime_error(
        "Can't calculate determinant, because of rows != columns");
  } else {
    S21Matrix triangle(*this);
    det = 1.0;
    for (int i = 0; i < rows_; i++) {
      int pivot = i;
      for (int j = i + 1; j < rows_; j++) {
        if (abs(triangle.matrix_[j][i]) > abs(triangle.matrix_[pivot][i])) {
          pivot = j;
        }
      }
      if (pivot != i) {
        triangle.swaprows(i, pivot);
        det *= -1;
      }
      if (triangle.matrix_[i][i] == 0) {
        det = 0;
      } else {
        det *= triangle.matrix_[i][i];
        for (int j = i + 1; j < rows_; j++) {
          double factor = triangle.matrix_[j][i] / triangle.matrix_[i][i];
          for (int k = i + 1; k < rows_; k++) {
            triangle.matrix_[j][k] -= factor * triangle.matrix_[i][k];
          }
        }
      }
    }
  }
  return det;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix buf(*this);
  if (rows_ != cols_)
    std::runtime_error(
        "Can't calculate Complemets Matrix, because of rows != columns");
  S21Matrix result(rows_, cols_);
  S21Matrix minor;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < rows_; j++) {
      result.matrix_[i][j] = MinorMatrixDet(i, j);
      result.matrix_[i][j] *= pow((-1.0), (i + j));
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = this->Determinant();

  if (det == 0) {
    std::runtime_error(
        "Can't calculate Inverse Matrix, because of Determinant = 0");
  }

  S21Matrix result;
  result = this->CalcComplements();
  result = result.Transpose();
  result *= (1.0 / det);
  return result;
}

// Operators

bool S21Matrix::operator==(const S21Matrix other) { return EqMatrix(other); }

S21Matrix& S21Matrix::operator=(const S21Matrix other) {
  if (*this == other) {
    ;
  }
  if ((rows_ == other.rows_) && (cols_ == other.cols_)) {
    CopyMatrix(other);
  } else {
    deallocation();
    rows_ = other.rows_;
    cols_ = other.cols_;
    InitializeMatrix();
    CopyMatrix(other);
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix buf(*this);
  buf.SumMatrix(other);
  return buf;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix buf(*this);
  buf.SubMatrix(other);
  return buf;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix buf(*this);
  buf.MulMatrix(other);
  return buf;
}

S21Matrix operator*(const double num, S21Matrix& other) { return other * num; }

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix buf(*this);
  buf.MulNumber(num);
  return buf;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return *this;
}

double& S21Matrix::operator()(int rows, int cols) {
  if (rows_ <= rows || cols_ <= cols || rows < 0 || cols < 0) {
    throw std::runtime_error("Wrong coordinates");
  }
  return matrix_[rows][cols];
}

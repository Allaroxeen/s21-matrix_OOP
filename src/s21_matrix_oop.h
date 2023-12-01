#ifndef S21_MATRIX_OOP_H
#define S21_MATRIX_OOP_H

#define Error_Rate 1e-6
#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  // Attributes

  int rows_, cols_;
  double** matrix_;

  // private methods

  void InitializeMatrix();
  void CopyMatrix(const S21Matrix& other);
  void deallocation();

  double MinorMatrixDet(const int a, const int b);

 public:
  void swaprows(int a, int b);
  // constructors

  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();  // Destructor

  // setters/getters

  void SetRows(int rows);
  void SetCols(int cols);
  int GetRows();
  int GetCols();

  // public methods

  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  double& operator()(int row, int col);
  bool operator==(const S21Matrix other);
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();

  // operators

  S21Matrix& operator=(const S21Matrix other);
  S21Matrix operator+(const S21Matrix&);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator*(const double num);
  friend S21Matrix operator*(const double num, const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  S21Matrix& operator*=(const S21Matrix& other);

  void PrintMatrix();
};

#endif  // S21MATRIX_OOP_H
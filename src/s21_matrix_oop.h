#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

class S21Matrix {
 public:
  int get_rows() const;
  int get_cols() const;

  void set_row(int rows);
  void set_col(int cols);
  S21Matrix(int rows, int cols);
  // конструктор по умолчанию
  // вызывается при создании объекта класса
  S21Matrix();
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  // operations
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix transpose() const;
  S21Matrix CalcComplements();
  double determinant();
  S21Matrix InverseMatrix();

  // перегрузка операторов
  bool operator==(const S21Matrix& other);
  bool operator!=(const S21Matrix& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double other);
  S21Matrix& operator+(const S21Matrix& other);
  S21Matrix& operator-(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  friend S21Matrix operator*(const S21Matrix& other, const double other1);
  friend S21Matrix operator*(const double other1, const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  double& operator()(int i, int j);
  double operator()(int i, int j) const;

  // functions for test
  S21Matrix(int rows, int cols, double data[]);
  void ToZero() const;

 private:
  int rows_, cols_;
  double** matrix_;

  // additional functions
  void PrintMatrix() const;
  void PrintConstructorU();
  S21Matrix GetMinor(int rows, int cols) const;
  void CreateMatrix();
  void Clear();
};
#endif  //  SRC_S21_MATRIX_OOP_H_

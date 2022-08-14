#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(constructors, T1) {
  S21Matrix one;
  S21Matrix second(2, 2);
  S21Matrix sec = std::move(second);
  S21Matrix third(sec);
  ASSERT_TRUE(third == sec);
}

TEST(equals_copy, T2) {
  double data1[4] = {1, 2, 3, 4};
  S21Matrix one(2, 2, data1);
  S21Matrix second;
  second = one;
  S21Matrix res(2, 2, data1);
  ASSERT_TRUE(one == res);
  ASSERT_TRUE(second == res);
}

TEST(matrix_sum, T1) {
  double data1[4] = {1, 2, 3, 4};
  double data2[4] = {1, 2, 2, 4};
  double data3[4] = {2, 4, 5, 8};
  S21Matrix one(2, 2, data1);
  S21Matrix second(2, 2, data2);
  S21Matrix res(2, 2, data3);
  one.SumMatrix(second);
  ASSERT_TRUE(one == res);
}

TEST(sum_throw, T1) {
  double data1[4] = {1, 2, 3, 4};
  double data2[9] = {1, 2, 2, 4, 5, 6, 7, 8, 9};
  double data3[6] = {1, 2, 3, 4, 5, 6};
  S21Matrix A(2, 2, data1);
  S21Matrix B(3, 3, data2);
  S21Matrix A1(3, 2, data3);
  EXPECT_THROW(A.SumMatrix(B), std::out_of_range);
  EXPECT_THROW(B.SubMatrix(A), std::out_of_range);
  EXPECT_THROW(B.MulMatrix(A), std::logic_error);
  EXPECT_THROW(A1.determinant(), std::logic_error);
  EXPECT_THROW(A1.CalcComplements(), std::logic_error);
}

TEST(matrix_sub, T1) {
  double data1[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  double data2[9] = {2, 2, 2, 2, 2, 2, 2, 2, 2};
  double data3[9] = {3, 4, 5, 6, 7, 8, 9, 10, 11};
  double data4[9] = {-1, 0, 1, 2, 3, 4, 5, 6, 7};
  S21Matrix one(3, 3, data1);
  S21Matrix second(3, 3, data2);
  S21Matrix third(3, 3, data3);
  S21Matrix fourth(3, 3, data4);
  S21Matrix res = third - second;
  ASSERT_TRUE(one == res);
  third.SubMatrix(second);
  ASSERT_TRUE(third == fourth);
  second.ToZero();
  res -= one;
  ASSERT_TRUE(res == second);
}

TEST(MulNumber, T1) {
  double data1[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  double data2[9] = {3, 6, 9, 12, 15, 18, 21, 24, 27};
  S21Matrix one(3, 3, data1);
  S21Matrix second(3, 3, data2);
  ASSERT_TRUE(second == one * 3);
  one.MulNumber(3);
  ASSERT_TRUE(second == one);
  one *= 1;
  ASSERT_TRUE(one == second);
}

TEST(MulMatrix, T1) {
  double data[] = {3, 2, 3, 4, 2, 7, 4, 4, 4};
  double data1[] = {1, 0, 2, 1, 1, 1, 4, 4, 4};
  double result[] = {17, 14, 20, 34, 30, 38, 24, 20, 28};
  S21Matrix one(3, 3, data);
  double a = one(0, 0);
  double b = 3;
  ASSERT_TRUE(a == b);
  S21Matrix second(3, 3, data1);
  S21Matrix res(3, 3, result);
  S21Matrix check(3, 3);
  check = one * second;
  ASSERT_TRUE(check == res);
  one.MulMatrix(second);
  ASSERT_TRUE(one == res);
  res.ToZero();
  one *= res;
  ASSERT_TRUE(one == res);
}

TEST(inverse_matrix_determinant_calc_complements_and_transpose, T1) {
  double data[] = {2, 4, 7, 14, 0, 5, 0, 1, 1};
  double data1[] = {-0.15625, 0.09375, 0.625,   -0.4375, 0.0625,
                    2.75,     0.4375,  -0.0625, -1.75};
  S21Matrix one(3, 3, data);
  S21Matrix second(3, 3, data1);
  one = one.InverseMatrix();
  ASSERT_TRUE(one == second);
}

TEST(inverse_and_equal, T1) {
  double data1[4] = {5, 0, 3.002, -4};
  double data2[4] = {5, 0, 3.003, -4};

  double data3[4] = {5, 0, 3, -4};
  double res[4] = {0.2, 0, 0.15, -0.25};

  S21Matrix a1(2, 2, data1);
  S21Matrix a2(2, 2, data2);
  S21Matrix a3(2, 2, data3);
  S21Matrix a4(2, 2, res);

  bool ravno = a1.EqMatrix(a2);
  ASSERT_EQ(ravno, 0);

  S21Matrix res1 = a3.InverseMatrix();
  ASSERT_TRUE(res1 == a4);
}

TEST(setters, T1) {
  double data1[4] = {3, 0, 0, 0};
  S21Matrix A(2, 2, data1);
  EXPECT_THROW(A.set_row(0), std::logic_error);
}

TEST(sets, T1) {
  double data[] = {3, 2, 3, 4, 2, 7, 4, 4, 4};
  double data1[] = {3, 2, 3, 4, 2, 7, 4, 4, 4, 0, 0, 0};
  double data2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  double data3[] = {1, 2, 3, 4};
  double data4[] = {1, 2, 0, 3, 4, 0};
  double data5[] = {3, 2, 3};
  S21Matrix one(3, 3, data);
  S21Matrix second(4, 3, data1);
  S21Matrix th(3, 3, data2);
  S21Matrix four(2, 2, data3);
  S21Matrix five(2, 3, data4);
  S21Matrix six(1, 3, data5);
  one.set_row(4);
  ASSERT_TRUE(one == second);
  th.set_col(3);
  ASSERT_TRUE(th == th);
  four.set_col(3);
  ASSERT_TRUE(four == five);
  int gr = th.get_rows();
  int b = 3;
  ASSERT_TRUE(gr == b);
  int gc = th.get_cols();
  ASSERT_TRUE(gc == b);
  one.set_row(1);
  ASSERT_TRUE(one == six);
}

int main() {
  testing::InitGoogleTest();
  return RUN_ALL_TESTS();
}

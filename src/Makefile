STD=-std=c++17
CFLAGS := -Wall -Werror -Wextra
LIBS := -lgtest
REPORTDIR := report
GCOV := --coverage 

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	FLAGS = -lgtest
	GCOV = -coverage
endif

ifeq ($(UNAME),Linux)
	FLAGS = -lgtest -lpthread -lrt -lm
	GCOV = -fprofile-arcs -ftest-coverage
endif

all: clean test leaks

s21_matrix_oop.a:
	g++ $(STD) -c S21_matrix_oop.cc -o s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o

test: s21_matrix_oop.a
	g++ $(CFLAGS) $(STD) $(GCOV) -c s21_matrix_oop.cc $(LIBS)
	g++ $(CFLAGS) $(STD) -c test.cc $(LIBS)
	g++ $(CFLAGS) $(STD) $(GCOV) -o test test.o s21_matrix_oop.o $(LIBS)
	./test

leaks: test
	leaks -atExit -- ./test

check:
	@cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	@python3 ../materials/linters/cpplint.py --extension=cpp *.cc
	@cppcheck *.cc
	@python3 ../materials/linters/cpplint.py --extension=h *.h
	@rm -rf CPPLINT.cfg

gcov_report:
	lcov -t "Unit-tests of s21_matrix_oop" -o s21_matrix_oop.info -c -d .
	genhtml -o $(REPORTDIR) s21_matrix_oop.info
	open -a "Google Chrome" $(REPORTDIR)/index.html
clean:
	@rm -rf *.o *.a *.gcna *.gcda *.gcno *.info report
	@rm -rf ./report
	@rm -rf test
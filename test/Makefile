CXX=g++
CXXFLAGS=-O0
CXXFLAGS+=-Wall -Wpedantic -Weffc++ -Warray-bounds -std=c++11 -g -pthread
CXXFLAGS+=-DPARALLEL_STRASSEN
CXXFLAGS+=-DTEST_MODE
LDFLAGS=-lgtest

TARGET=qmatrix_test

INCLUDE=-I../

all: ${TARGET}

${TARGET}: matrix_strassen.o MatrixTest.o main.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} matrix_strassen.o MatrixTest.o main.o -o ${TARGET}

matrix_strassen.o: ../matrix_strassen.cpp
	${CXX} ${INCLUDE} ${CXXFLAGS} -c ../matrix_strassen.cpp

MatrixTest.o: MatrixTest.cpp
	${CXX} ${INCLUDE} ${CXXFLAGS} -c MatrixTest.cpp

main.o: main.cpp
	${CXX} ${INCLUDE} ${CXXFLAGS} -c main.cpp

clean:
	rm -rf *.o ${TARGET}

CXX=g++
CXXFLAGS=-O3
CXXFLAGS+=-Wall -Wpedantic -Weffc++ -Warray-bounds -std=c++11 -g -pthread
#CXXFLAGS+=-DTRIVIAL_ALGORITHM
CXXFLAGS+=-DPARALLEL_STRASSEN
#CXXFLAGS+=-DTEST_MODE
CXXFLAGS+=-ftree-vectorize -msse2 -ftree-vectorizer-verbose=5

TARGET=qmatrix

all: ${TARGET}

${TARGET}: matrix_strassen.o main.o
	${CXX} ${CXXFLAGS} matrix_strassen.o main.o -o ${TARGET}

matrix_strassen.o: matrix_strassen.cpp
	${CXX} ${CXXFLAGS} -c matrix_strassen.cpp

main.o: main.cpp
	${CXX} ${CXXFLAGS} -c main.cpp

clean:
	rm -rf *.o ${TARGET}

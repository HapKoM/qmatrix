CXX=g++
CXXFLAGS=-O0
CXXFLAGS+=-Wall -Wpedantic -Weffc++ -Warray-bounds -std=c++11 -g -pthread
CXXFLAGS+=-DPARALLEL_STRASSEN

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

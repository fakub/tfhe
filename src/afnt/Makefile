

all:
	g++ -Wall -c -Iinclude main.cpp
	g++ -Wall nayuki/portable/fft-x8664-avx-aux.c.o nayuki/portable/fft-model-of-x8664-avx.c.o main.o -o main

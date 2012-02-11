all: generate

generate: generate.c
	gcc -Wall -std=c99 -o generate -lm -ltiff -fopenmp -O2 generate.c

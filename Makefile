LDFLAGS += -lm -ltiff
CFLAGS += -Wall -Werror -std=c99 -fopenmp -O3
all: generate

generate: generate.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

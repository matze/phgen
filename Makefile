LDFLAGS += -lm -ltiff
CFLAGS += -Wall -Werror -std=c99 -fopenmp -O3
BIN = generate

.PHONY: clean

all: $(BIN)

$(BIN): generate.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(BIN)

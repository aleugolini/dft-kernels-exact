CC ?= cc
# -std=gnu99: older gcc versions (e.g. 4.8) default to C89. Add -march=native for speed
# when the program is built on the machine that runs it.
CFLAGS ?= -O3 -std=gnu99 -Wall

fs_search: fs_search.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f fs_search

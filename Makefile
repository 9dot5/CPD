CC = gcc
CFLAGS = -Wall -O2

life3d: life3d.c
	$(CC) $(CFLAGS) -o $@ life3d.c

clean:
	rm -f life3d

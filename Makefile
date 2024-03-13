CC = gcc
CFLAGS = -O2 -fopenmp

life3d: life3d.c
	$(CC) $(CFLAGS) -o life3d life3d.c

clean:
	rm -f life3d
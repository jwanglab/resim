CC     = gcc
LIBS   = -lz -lm

OBJECTS = resim

all: $(OBJECTS)

resim:
	$(CC) $(CFLAGS) -o resim main.c cmap.c bnx.c sim.c digest.c $(LIBS)

.PHONY: clean
clean:
	-rm $(OBJECTS)

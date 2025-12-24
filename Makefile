CC = gcc

# macOS with Homebrew LAPACK
CFLAGS = -O2 -Wall -std=c11 -I/opt/homebrew/include -I/opt/homebrew/opt/lapack/include
LIBS = -L/opt/homebrew/opt/lapack/lib/ -llapacke -llapack -lblas -lm

SRC = main.c gendata.c ula.c music.c
OBJ = $(SRC:.c=.o)
TARGET = music

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJ) $(TARGET)

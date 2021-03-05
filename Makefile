# Makefile example
# Authors:
#	Alessandro Budroni, May 2019

# directories variables
OBJ_dir=build
SRC_dir=src
INCLUDE_dir=include

# global variables
CC=gcc
MD=mkdir
CFLAGS= -std=c11 -O3 -g #-m64 -Wformat=0 -Wno-unused-function -Wno-unused-result -D_FILE_OFFSET_BITS=64 -DDEBUG -D_DEBUG -pedantic
IDIR = -I /usr/include -I /usr/local/include/ -I ./$(INCLUDE_dir)
LDIR = -L /usr/lib/ -L /usr/local/lib/ -L ./$(OBJ_dir) 
LIBS= -lm

# headers
_HEADER_files = config.h  utils.h
HEADER_files = $(patsubst %,$(INCLUDE_dir)/%,$(_HEADER_files))
# source and objects
_SRC_files =  utils.c
SRC_files = $(patsubst %,$(SRC_dir)/%,$(_SRC_files))
OBJ_files = $(patsubst %.c,$(OBJ_dir)/%.o,$(_SRC_files))

all: $(OBJ_dir) main

# create build directory
$(OBJ_dir):
	$(MD) $@

# compile
$(OBJ_dir)/%.o: $(SRC_dir)/%.c
	$(CC) $(CFLAGS) $< -c $(IDIR) $(LDIR) -o $@ $(LIBS)

main: $(OBJ_files)
	$(CC) $(CFLAGS) $(SRC_dir)/$@.c -o $(OBJ_dir)/$@ $(IDIR) $(LDIR) $^ $(LIBS) 

.PHONY: clean all

clean:
	rm -rf $(OBJ_dir)

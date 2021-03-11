# Makefile example
# Authors:
#	Alessandro Budroni, May 2019

# directories variables
OBJ_dir=build
SRC_dir=src
TEST_dir=test
INCLUDE_dir=include

# global variables
CC=gcc
MD=mkdir
CFLAGS= -std=c11 -O3 # -fsanitize=address -g -fno-omit-frame-pointer # -m64 -Wformat=0 -Wno-unused-function -Wno-unused-result -D_FILE_OFFSET_BITS=64 -DDEBUG -D_DEBUG -pedantic
IDIR = -I /usr/include -I /usr/local/include/ -I ./$(INCLUDE_dir)
LDIR = -L /usr/lib/ -L /usr/local/lib/ -L ./$(OBJ_dir) 
LIBS= -lm

# headers
_HEADER_files = config.h  utils.h  lwe_instance.h transition_times2_modq.h position_values_2_category_index.h transition_bkw_step_smooth_lms.h lookup_tables.h transition_bkw_step_final.h solve_fwht.h
HEADER_files = $(patsubst %,$(INCLUDE_dir)/%,$(_HEADER_files))
# source and objects
_SRC_files =  utils.c lwe_instance.c transition_times2_modq.c position_values_2_category_index.c transition_bkw_step_smooth_lms.c lookup_tables.c transition_bkw_step_final.c solve_fwht.c
SRC_files = $(patsubst %,$(SRC_dir)/%,$(_SRC_files))
OBJ_files = $(patsubst %.c,$(OBJ_dir)/%.o,$(_SRC_files))

all: $(OBJ_dir) test40_005 test20_005 test20_005_large

# create build directory
$(OBJ_dir):
	$(MD) $@

# compile
$(OBJ_dir)/%.o: $(SRC_dir)/%.c
	$(CC) $(CFLAGS) $< -c $(IDIR) $(LDIR) -o $@ $(LIBS)

test40_005: $(OBJ_files)
	$(CC) $(CFLAGS) $(TEST_dir)/$@.c -o $(OBJ_dir)/$@ $(IDIR) $(LDIR) $^ $(LIBS)

test20_005: $(OBJ_files)
	$(CC) $(CFLAGS) $(TEST_dir)/$@.c -o $(OBJ_dir)/$@ $(IDIR) $(LDIR) $^ $(LIBS)

test20_005_large: $(OBJ_files)
	$(CC) $(CFLAGS) $(TEST_dir)/$@.c -o $(OBJ_dir)/$@ $(IDIR) $(LDIR) $^ $(LIBS)

.PHONY: clean all

clean:
	rm -rf $(OBJ_dir)

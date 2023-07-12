###############################################################################
# Makefile
#
# AUTHOR  : Felipe Belem
# DATE    : 2021-02-24
# LICENSE : MIT License
# EMAIL   : felipe.belem@ic.unicamp.br
###############################################################################
#==============================================================================
# COMPILING SETTINGS
#==============================================================================
CC = gcc #mingw-64 for Windows, or clang for MacOS
CFLAGS= -std=gnu11 -fPIC

# Directories -----------------------------------------------------------------
HOME_DIR = .
INC_DIR = $(HOME_DIR)/include
SRC_DIR = $(HOME_DIR)/src
OBJ_DIR = $(HOME_DIR)/obj
LIB_DIR = $(HOME_DIR)/lib
BIN_DIR = $(HOME_DIR)/bin
DEMO_DIR = $(HOME_DIR)/demo

INCS = -I$(INC_DIR)
LIBS_LD = -L$(LIB_DIR)
LIB_NAME = phd
LIBS_LINK = -l$(LIB_NAME)
LIB_SFIX = .a#.lib for Windows

# Compiler --------------------------------------------------------------------
# Enabling debugging GNU GCC flags
IFT_DEBUG = NO

# Enabling external libraries
# OpenMP: Requires version 4.5
# LibPNG: Requires version 1.6.29
# LibJPEG: Requires version 9.0
IFT_OMP = YES
IFT_LIBPNG = YES
IFT_LIBJPEG = YES

# It is expecting a GNU GCC compiler. For other compilers, modifications 
# might be necessary
ifeq ($(IFT_DEBUG),YES)
	CFLAGS += -Og -g -pedantic -ggdb -pg -Wfatal-errors -Wall -Wextra \
						-Wno-unused-parameter -DIFT_DEBUG
else
	CFLAGS += -O3

	ifeq ($(IFT_OMP), YES)
		CFLAGS += -fopenmp -DIFT_OMP
		LIBS_LINK += -lgomp
	endif
endif

ifeq ($(IFT_LIBJPEG),YES)
	# If you desire to indicate another library version to be used (whether it
	# is a shared or static one), update the following 2 variables accordingly
	INCS += -I$(HOME_DIR)/externals/libjpeg/include# /path/to/libjpeg/headers
	LIBS_LD += -L$(HOME_DIR)/externals/libjpeg/lib# /path/to/libjpeg/library
	LIBS_LINK += -ljpeg
	CFLAGS += -DIFT_LIBJPEG
endif

ifeq ($(IFT_LIBPNG),YES) 
	# If you desire to indicate another library version to be used (whether it
	# is a shared or static one), update the following 2 variables accordingly
	INCS += -I$(HOME_DIR)/externals/libpng/include# /path/to/libpng/headers
	LIBS_LD += -L$(HOME_DIR)/externals/libpng/lib# /path/to/libpng/library
	LIBS_LINK += -lpng
	CFLAGS += -DIFT_LIBPNG

	# ZLib: Requires version 1.2.11
	# If you desire to indicate another library version to be used (whether it
 	# is a shared or static one), update the following 2 variables accordingly
	INCS += -I$(HOME_DIR)/externals/libpng/zlib/include# /path/to/zlib/headers
	LIBS_LD += -L$(HOME_DIR)/externals/libpng/zlib/lib# /path/to/zlib/library
	LIBS_LINK += -lz
endif

LIBS_LINK += -lm

# Files -----------------------------------------------------------------------
SRC_FILES = $(wildcard $(SRC_DIR)/*.c)
DEMO_FILES = $(wildcard $(DEMO_DIR)/*.c)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRC_FILES))

#==============================================================================
# RULES
#==============================================================================
.PHONY: all clean demo lib obj refresh tidy

all: lib demo

# Compiling -------------------------------------------------------------------
status:
	@echo "\nCompilation Flags ------------------------------------------"
	@echo "- IFT_DEBUG: $(IFT_DEBUG)"
	@echo "- IFT_OMP: $(IFT_OMP)"
	@echo "- IFT_LIBPNG: $(IFT_LIBPNG)"
	@echo "- IFT_LIBJPEG: $(IFT_LIBJPEG)"
	@echo "------------------------------------------------------------\n"

lib: status obj
	@mkdir -p $(LIB_DIR)
	ar csr $(LIB_DIR)/lib$(LIB_NAME)$(LIB_SFIX) $(OBJ_FILES)
	@echo "\n----- Library was successfully built\n"

obj:
	@make -j $(OBJ_FILES)
	@echo "\n----- All sources were compiled\n"

demo:
	@make -j $(patsubst %.c, %, $(DEMO_FILES))
	@echo "\n----- All demos were compiled\n"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@ $(LIBS_LD) $(LIBS_LINK)

$(DEMO_DIR)/%:
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $(INCS) $@.c -o $(BIN_DIR)/$(@F) $(LIBS_LD) $(LIBS_LINK)

# Cleaning --------------------------------------------------------------------
tidy:
	$(RM) -r $(BIN_DIR)

clean: tidy
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(LIB_DIR)

refresh: clean lib

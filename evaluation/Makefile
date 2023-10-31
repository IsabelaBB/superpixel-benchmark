#==============================================================================
# Paths
#==============================================================================
DEMO_DIR = .

BIN_DIR = $(DEMO_DIR)/bin
STB_DIR = $(DEMO_DIR)/externals
SRC_DIR = $(DEMO_DIR)/src
LIB_DIR = $(DEMO_DIR)/lib
INCLUDE_DIR = $(DEMO_DIR)/include
OBJ_DIR = $(DEMO_DIR)/obj

CC = gcc -g
CFLAGS = -Wall -fPIC -std=gnu11 -pedantic -Wno-unused-result -O3 -fopenmp 

CXX = g++ -g
DEMOFLAGS = -Wall -fPIC -pedantic -Wno-unused-result -O3 -fopenmp
CXXFLAGS = -O3 -lstdc++ -fPIC -ffast-math -march=skylake -mfma -Wall -Wno-unused-result

LIBS = -lm `pkg-config --cflags --libs opencv4`

HEADER_INC = -I $(STB_DIR) -I $(INCLUDE_DIR)

LIBS_LD = -L$(LIB_DIR)
LIB_NAME = supColorVaricance
LIBS_LINK = -l$(LIB_NAME)
LIB_SFIX = .a#.lib for Windows

IFT_LIBPNG = YES
IFT_LIBJPEG = YES

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

#==============================================================================
# Rules
#==============================================================================
.PHONY: all c clean lib

all: folders lib c 

folders: 
	@mkdir -p $(OBJ_DIR) 
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(LIB_DIR)

lib: obj
	$(eval ALL_OBJS := $(wildcard $(OBJ_DIR)/*.o))
	ar csr $(LIB_DIR)/lib$(LIB_NAME)$(LIB_SFIX) $(ALL_OBJS)

obj: \
	$(OBJ_DIR)/Utils.o \
	$(OBJ_DIR)/Color.o \
	$(OBJ_DIR)/PrioQueue.o \
	$(OBJ_DIR)/Image.o \
	$(OBJ_DIR)/Eval.o \
	$(OBJ_DIR)/ift.o 
	

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INCLUDE_DIR)/%.h
	$(CC) $(CFLAGS) -c $< -o $@ $(HEADER_INC) $(LIBS_LD) $(LIBS_LINK) $(LIBS)

c: lib
	$(CXX) $(DEMOFLAGS) $(CXXFLAGS) ./main.cpp -o $(BIN_DIR)/main $(HEADER_INC) $(LIBS_LD) $(LIBS_LINK) $(LIBS)

clean:
	rm -rf $(OBJ_DIR)/ ;
	rm -rf $(BIN_DIR)/ ;
	rm -rf $(LIB_DIR)/ ;

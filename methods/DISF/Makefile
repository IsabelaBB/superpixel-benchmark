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

CC = gcc
CFLAGS = -Wall -fPIC -std=gnu11 -pedantic -Wno-unused-result -O3 -fopenmp
LIBS = -lm

HEADER_INC = -I $(STB_DIR) -I $(INCLUDE_DIR)
LIB_INC = -L $(LIB_DIR) -ldisf

#==============================================================================
# Rules
#==============================================================================
.PHONY: all c clean lib folders

all: folders lib c

folders: 
	@mkdir -p $(OBJ_DIR) 
	@mkdir -p $(BIN_DIR)

lib: obj
	$(eval ALL_OBJS := $(wildcard $(OBJ_DIR)/*.o))
	@mkdir -p $(LIB_DIR)
	ar csr $(LIB_DIR)/libdisf.a $(ALL_OBJS) ;

obj: \
	$(OBJ_DIR)/Utils.o \
	$(OBJ_DIR)/IntList.o \
	$(OBJ_DIR)/Color.o \
	$(OBJ_DIR)/PrioQueue.o \
	$(OBJ_DIR)/Image.o \
	$(OBJ_DIR)/DISF.o 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INCLUDE_DIR)/%.h
	@mkdir -p $(@D) 
	$(CC) $(CFLAGS) -c $< -o $@ $(HEADER_INC) $(LIBS)

c: lib
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) DISF_demo.c -o $(BIN_DIR)/DISF_demo $(HEADER_INC) $(LIB_INC) $(LIBS)

clean:
	rm -rf $(OBJ_DIR) ;
	rm -rf $(LIB_DIR) ;
	rm -rf $(BIN_DIR) ;

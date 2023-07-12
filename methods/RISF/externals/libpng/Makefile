# makefile for libpng using gcc (generic, static library)

# This code is released under the libpng license.
# For conditions of distribution and use, see the disclaimer
# and license in png.h


# Location of the zlib library and include files
ZLIBINC = ../zlib/include
ZLIBLIB = ../zlib/lib



LIB=./lib
INCLUDESDIR = ./include
SRC=./src
OBJ=./obj
BIN=./bin
#CC= cc

INCLUDE = -I$(INCLUDESDIR)

WARNMORE = -Wwrite-strings -Wpointer-arith -Wshadow \
	-Wmissing-declarations -Wtraditional -Wcast-align \
	-Wstrict-prototypes -Wmissing-prototypes # -Wconversion
CPPFLAGS = -I$(ZLIBINC) # -DPNG_DEBUG=5
FLAGS = $(CPPFLAGS) $(CFLAGS)
LIBS = -lz -lm

WARNMORE = -Wwrite-strings -Wpointer-arith -Wshadow \
	-Wmissing-declarations -Wtraditional -Wcast-align \
	-Wstrict-prototypes -Wmissing-prototypes # -Wconversion
CPPFLAGS = -I$(ZLIBINC) # -DPNG_DEBUG=5
CFLAGS = -fPIC -I$(ZLIBINC) -W -Wall -std=gnu11 -O2 # $(WARNMORE) 
LDFLAGS =
LIBS = -lz -lm

# Pre-built configuration
# See scripts/pnglibconf.mak for more options
#PNGLIBCONF_H_PREBUILT = pnglibconf.h.prebuilt


libpng: zlib $(LIB)/libiftpng.a
	echo "libiftpng.a built..."

$(LIB)/libiftpng.a: \
$(OBJ)/png.o \
$(OBJ)/pngerror.o \
$(OBJ)/pngget.o \
$(OBJ)/pngmem.o \
$(OBJ)/pngpread.o \
$(OBJ)/pngread.o \
$(OBJ)/pngrio.o \
$(OBJ)/pngrtran.o \
$(OBJ)/pngrutil.o \
$(OBJ)/pngset.o \
$(OBJ)/pngtrans.o \
$(OBJ)/pngwio.o \
$(OBJ)/pngwrite.o \
$(OBJ)/pngwtran.o \
$(OBJ)/pngwutil.o \

	ar csr $(LIB)/libiftpng.a \
$(OBJ)/*.o

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) -c $< $(INCLUDE) -o $@

zlib:
	cd ../zlib; $(MAKE) ; cd -\

clean:
	rm -f $(LIB)/*.a; rm -f $(OBJ)/*.o;

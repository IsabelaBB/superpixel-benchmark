# The name of your C compiler:
#CC= cc
CFLAGS= -fPIC -O3 -std=gnu11
LIB=./lib
INCLUDESDIR  = ./include
SRC=./src
OBJ=./obj
BIN=./bin

INCLUDE = -I$(INCLUDESDIR)

# Put here the object file name for the correct system-dependent memory
# manager file.  For Unix this is usually jmemnobs.o, but you may want
# to use jmemansi.o or jmemname.o if you have limited swap space.
SYSDEPMEM= jmemnobs

libjpeg: $(LIB)/libjpeg.a
	echo "libjpeg.a built..."

$(LIB)/libjpeg.a: \
$(OBJ)/jcapimin.o \
$(OBJ)/jcapistd.o \
$(OBJ)/jcarith.o \
$(OBJ)/jctrans.o \
$(OBJ)/jcparam.o \
$(OBJ)/jdatadst.o \
$(OBJ)/jcinit.o \
$(OBJ)/jcmaster.o \
$(OBJ)/jcmarker.o \
$(OBJ)/jcmainct.o \
$(OBJ)/jcprepct.o \
$(OBJ)/jccoefct.o \
$(OBJ)/jccolor.o \
$(OBJ)/jcsample.o \
$(OBJ)/jchuff.o \
$(OBJ)/jcdctmgr.o \
$(OBJ)/jfdctfst.o \
$(OBJ)/jfdctflt.o \
$(OBJ)/jfdctint.o \
$(OBJ)/jdapimin.o \
$(OBJ)/jdapistd.o \
$(OBJ)/jdarith.o \
$(OBJ)/jdtrans.o \
$(OBJ)/jdatasrc.o \
$(OBJ)/jdmaster.o \
$(OBJ)/jdinput.o \
$(OBJ)/jdmarker.o \
$(OBJ)/jdhuff.o \
$(OBJ)/jdmainct.o \
$(OBJ)/jdcoefct.o \
$(OBJ)/jdpostct.o \
$(OBJ)/jddctmgr.o \
$(OBJ)/jidctfst.o \
$(OBJ)/jidctflt.o \
$(OBJ)/jidctint.o \
$(OBJ)/jdsample.o \
$(OBJ)/jdcolor.o \
$(OBJ)/jquant1.o \
$(OBJ)/jquant2.o \
$(OBJ)/jdmerge.o \
$(OBJ)/jaricom.o \
$(OBJ)/jcomapi.o \
$(OBJ)/jutils.o \
$(OBJ)/jerror.o \
$(OBJ)/jmemmgr.o \
$(OBJ)/$(SYSDEPMEM).o \


	ar csr $(LIB)/libjpeg.a \
$(OBJ)/*.o

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) -c $< $(INCLUDE) -o $@

#compression library object files:
#jcapimin.c jcapistd.c jcarith.c jctrans.c jcparam.c
#jdatadst.c jcinit.c jcmaster.c jcmarker.c jcmainct.c
#jcprepct.c jccoefct.c jccolor.c jcsample.c jchuff.c
#jcdctmgr.c jfdctfst.c jfdctflt.c jfdctint.c

#decompression library object files:
#jdapimin.c jdapistd.c jdarith.c jdtrans.c jdatasrc.c
#jdmaster.c jdinput.c jdmarker.c jdhuff.c jdmainct.c
#jdcoefct.c jdpostct.c jddctmgr.c jidctfst.c jidctflt.c
#jidctint.c jdsample.c jdcolor.c jquant1.c jquant2.c
#jdmerge.c

# library object files common to compression and decompression:
#jaricom.c jcomapi.c jutils.c jerror.c jmemmgr.c

#System dependence memory:

#This file provides a really simple implementation of the system-
#dependent portion of the JPEG memory manager.  This implementation
#assumes that no backing-store files are needed: all required space
#can be obtained from malloc().
#This is very portable in the sense that it'll compile on almost anything,
#but you'd better have lots of main memory (or virtual memory) if you want
#to process big images.
#Note that the max_memory_to_use option is ignored by this implementation.
#jmemnobs.c

clean:
	rm -f $(LIB)/*.a; rm -f $(OBJ)/*.o;

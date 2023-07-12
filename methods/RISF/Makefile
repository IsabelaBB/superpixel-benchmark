CC = gcc
IDIR := include
SDIR := src
ODIR := obj
BIN := bin
DEMO := demo

# external libs
EXTERNALS_DIR = ./externals
LIBJPEG = $(EXTERNALS_DIR)/libjpeg/include
LIBPNG = $(EXTERNALS_DIR)/libpng
LIBZ = $(EXTERNALS_DIR)/zlib

INCLUDES = -I$(IDIR) -I$(LIBJPEG)/include -I$(LIBPNG)/include -I$(LIBZ)/include
LIBS = -L$(LIBJPEG)/lib -L$(LIBPNG)/lib -L$(LIBZ)/lib -ljpeg -liftpng -lm -liftz -lz
FLAGS= -O3 -DNDEBUG -Wall -Wno-unused-result -fopenmp -pthread -std=gnu11 -fPIC -pedantic $(INCLUDES)

#shared includes
_DEPS = iftSuperpixelFeatures.h \
	iftExperimentUtility.h \
	iftRISF.h \
	ift.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

#shared objects
_OBJ = iftSuperpixelFeatures.o \
       iftExperimentUtility.o \
       iftRISF.o \
       ift.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: external_libs $(BIN)/iftForceLabelMapConnectivity $(BIN)/iftMidLevelRegionMerge $(BIN)/iftOverlayBorders $(BIN)/iftRISF_segmentation $(BIN)/iftShowConnectivityProblems 
  
########### EXTERNAL LIBS ###########
$(LIBJPEG)/lib/libjpeg.a:
	cd $(LIBJPEG); $(MAKE) ; cd -\

$(LIBPNG)/lib/libiftpng.a:
	cd $(LIBPNG); $(MAKE) ; cd -\

$(LIBZ)/lib/libiftz.a:
	cd $(LIBZ); $(MAKE) ; cd -\

external_libs: $(LIBJPEG)/lib/libjpeg.a $(LIBPNG)/lib/libiftpng.a $(LIBZ)/lib/libiftz.a


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) $(FLAGS) -c -o $@ $< $(LIBS)


$(BIN)/iftComputeSegmentationMetrics: $(DEMO)/iftComputeSegmentationMetrics.c $(OBJ)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

$(BIN)/iftForceLabelMapConnectivity: $(DEMO)/iftForceLabelMapConnectivity.c $(OBJ)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

$(BIN)/iftMidLevelRegionMerge: $(DEMO)/iftMidLevelRegionMerge.c $(OBJ)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

$(BIN)/iftOverlayBorders: $(DEMO)/iftOverlayBorders.c $(OBJ)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

$(BIN)/iftRISF_segmentation: $(DEMO)/iftRISF_segmentation.c $(OBJ)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

$(BIN)/iftShowConnectivityProblems: $(DEMO)/iftShowConnectivityProblems.c $(OBJ)
	$(CC) $(FLAGS) $^ -o $@ $(LIBS)

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o
	rm -df $(ODIR)
	rm -f $(BIN)/*;

$(OBJ): $(ODIR)

$(ODIR):
	mkdir $(ODIR)

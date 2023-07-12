LIB=./lib
INCLUDE=./include
BIN=./bin
SRC=./src
OBJ=./obj
DEMO=./demo

#Compiladores
CXX=g++

#-march=native
FLAGS= -Wall -O3 -msse -fopenmp
LINKS= -lpthread -lz -lm 

all: folders libgft $(BIN)/superpixels $(BIN)/supervoxels
	@echo
	@echo "libgft.a built..."
	@echo "DONE."

folders: 
	@mkdir -p $(OBJ) 
	@mkdir -p $(BIN)
	@mkdir -p $(LIB)

libgft: $(LIB)/libgft.a

GFTLIB  = -L$(LIB) -lgft
GFTFLAGS  = -I$(INCLUDE)


$(LIB)/libgft.a: \
$(OBJ)/gft_common.o \
$(OBJ)/gft_adjrel.o \
$(OBJ)/gft_image32.o \
$(OBJ)/gft_image32f.o \
$(OBJ)/gft_image64.o \
$(OBJ)/gft_cimage.o \
$(OBJ)/gft_cimage32f.o \
$(OBJ)/gft_morphology.o \
$(OBJ)/gft_queue.o \
$(OBJ)/gft_gpqueue_by_Falcao.o \
$(OBJ)/gft_graph.o \
$(OBJ)/gft_sparsegraph.o \
$(OBJ)/gft_layeredgraph.o \
$(OBJ)/gft_heap.o \
$(OBJ)/gft_heap64f.o \
$(OBJ)/gft_heap_lex.o \
$(OBJ)/gft_ift.o \
$(OBJ)/gft_string.o \
$(OBJ)/gft_arraylist.o \
$(OBJ)/gft_filelist.o \
$(OBJ)/gft_color.o \
$(OBJ)/gft_highlight.o \
$(OBJ)/gft_superpixels.o \
$(OBJ)/gft_evaluation.o \
$(OBJ)/gft_scene.o \
$(OBJ)/gft_scene32.o \
$(OBJ)/gft_scene64.o \
$(OBJ)/gft_scene16.o \
$(OBJ)/gft_scene8.o \
$(OBJ)/gft_adjrel3.o \
$(OBJ)/gft_scnmath.o \
$(OBJ)/gft_matrix.o \
$(OBJ)/gft_bmap.o \
$(OBJ)/gft_set.o \
$(OBJ)/gft_stack.o \
$(OBJ)/gft_analysis.o \
$(OBJ)/gft_pqueue32.o \
$(OBJ)/gft_pqueue16.o \
$(OBJ)/gft_hashtable.o \
$(OBJ)/gft_curve.o \
$(OBJ)/gft_radiometric.o \
$(OBJ)/gft_radiometric3.o \
$(OBJ)/gft_filtering.o \
$(OBJ)/gft_filtering3.o \
$(OBJ)/gft_adjregion3.o \
$(OBJ)/gft_seedmap3.o \
$(OBJ)/gft_segmentation3.o \
$(OBJ)/gft_marker.o \
$(OBJ)/gft_features.o \
$(OBJ)/gft_bzlib.o \
$(OBJ)/nifti1_io.o \
$(OBJ)/znzlib.o
	ar csr $(LIB)/libgft.a \
$(OBJ)/gft_common.o \
$(OBJ)/gft_adjrel.o \
$(OBJ)/gft_image32.o \
$(OBJ)/gft_image32f.o \
$(OBJ)/gft_image64.o \
$(OBJ)/gft_cimage.o \
$(OBJ)/gft_cimage32f.o \
$(OBJ)/gft_morphology.o \
$(OBJ)/gft_queue.o \
$(OBJ)/gft_gpqueue_by_Falcao.o \
$(OBJ)/gft_graph.o \
$(OBJ)/gft_sparsegraph.o \
$(OBJ)/gft_layeredgraph.o \
$(OBJ)/gft_heap.o \
$(OBJ)/gft_heap64f.o \
$(OBJ)/gft_heap_lex.o \
$(OBJ)/gft_ift.o \
$(OBJ)/gft_string.o \
$(OBJ)/gft_arraylist.o \
$(OBJ)/gft_filelist.o \
$(OBJ)/gft_color.o \
$(OBJ)/gft_highlight.o \
$(OBJ)/gft_superpixels.o \
$(OBJ)/gft_evaluation.o \
$(OBJ)/gft_scene.o \
$(OBJ)/gft_scene32.o \
$(OBJ)/gft_scene64.o \
$(OBJ)/gft_scene16.o \
$(OBJ)/gft_scene8.o \
$(OBJ)/gft_adjrel3.o \
$(OBJ)/gft_scnmath.o \
$(OBJ)/gft_matrix.o \
$(OBJ)/gft_bmap.o \
$(OBJ)/gft_set.o \
$(OBJ)/gft_stack.o \
$(OBJ)/gft_analysis.o \
$(OBJ)/gft_pqueue32.o \
$(OBJ)/gft_pqueue16.o \
$(OBJ)/gft_hashtable.o \
$(OBJ)/gft_curve.o \
$(OBJ)/gft_radiometric.o \
$(OBJ)/gft_radiometric3.o \
$(OBJ)/gft_filtering.o \
$(OBJ)/gft_filtering3.o \
$(OBJ)/gft_adjregion3.o \
$(OBJ)/gft_seedmap3.o \
$(OBJ)/gft_segmentation3.o \
$(OBJ)/gft_marker.o \
$(OBJ)/gft_features.o \
$(OBJ)/gft_bzlib.o \
$(OBJ)/nifti1_io.o \
$(OBJ)/znzlib.o


$(OBJ)/gft_common.o: $(SRC)/gft_common.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_common.cpp -I$(INCLUDE) -o $(OBJ)/gft_common.o

$(OBJ)/gft_adjrel.o: $(SRC)/gft_adjrel.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_adjrel.cpp -I$(INCLUDE) -o $(OBJ)/gft_adjrel.o

$(OBJ)/gft_image32.o: $(SRC)/gft_image32.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_image32.cpp -I$(INCLUDE) -o $(OBJ)/gft_image32.o

$(OBJ)/gft_image32f.o: $(SRC)/gft_image32f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_image32f.cpp -I$(INCLUDE) -o $(OBJ)/gft_image32f.o

$(OBJ)/gft_image64.o: $(SRC)/gft_image64.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_image64.cpp -I$(INCLUDE) -o $(OBJ)/gft_image64.o

$(OBJ)/gft_cimage.o: $(SRC)/gft_cimage.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_cimage.cpp -I$(INCLUDE) -o $(OBJ)/gft_cimage.o

$(OBJ)/gft_cimage32f.o: $(SRC)/gft_cimage32f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_cimage32f.cpp -I$(INCLUDE) -o $(OBJ)/gft_cimage32f.o

$(OBJ)/gft_morphology.o: $(SRC)/gft_morphology.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_morphology.cpp -I$(INCLUDE) -o $(OBJ)/gft_morphology.o

$(OBJ)/gft_queue.o: $(SRC)/gft_queue.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_queue.cpp -I$(INCLUDE) -o $(OBJ)/gft_queue.o

$(OBJ)/gft_gpqueue_by_Falcao.o: $(SRC)/gft_gpqueue_by_Falcao.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_gpqueue_by_Falcao.cpp -I$(INCLUDE) -o $(OBJ)/gft_gpqueue_by_Falcao.o

$(OBJ)/gft_graph.o: $(SRC)/gft_graph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_graph.cpp -I$(INCLUDE) -o $(OBJ)/gft_graph.o

$(OBJ)/gft_sparsegraph.o: $(SRC)/gft_sparsegraph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_sparsegraph.cpp -I$(INCLUDE) -o $(OBJ)/gft_sparsegraph.o

$(OBJ)/gft_layeredgraph.o: $(SRC)/gft_layeredgraph.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_layeredgraph.cpp -I$(INCLUDE) -o $(OBJ)/gft_layeredgraph.o

$(OBJ)/gft_heap.o: $(SRC)/gft_heap.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap.o

$(OBJ)/gft_heap64f.o: $(SRC)/gft_heap64f.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap64f.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap64f.o

$(OBJ)/gft_heap_lex.o: $(SRC)/gft_heap_lex.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_heap_lex.cpp -I$(INCLUDE) -o $(OBJ)/gft_heap_lex.o

$(OBJ)/gft_ift.o: $(SRC)/gft_ift.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_ift.cpp -I$(INCLUDE) -o $(OBJ)/gft_ift.o

$(OBJ)/gft_string.o: $(SRC)/gft_string.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_string.cpp -I$(INCLUDE) -o $(OBJ)/gft_string.o

$(OBJ)/gft_arraylist.o: $(SRC)/gft_arraylist.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_arraylist.cpp -I$(INCLUDE) -o $(OBJ)/gft_arraylist.o

$(OBJ)/gft_filelist.o: $(SRC)/gft_filelist.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_filelist.cpp -I$(INCLUDE) -o $(OBJ)/gft_filelist.o

$(OBJ)/gft_highlight.o: $(SRC)/gft_highlight.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_highlight.cpp -I$(INCLUDE) -o $(OBJ)/gft_highlight.o

$(OBJ)/gft_color.o: $(SRC)/gft_color.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_color.cpp -I$(INCLUDE) -o $(OBJ)/gft_color.o


$(OBJ)/gft_superpixels.o: $(SRC)/gft_superpixels.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_superpixels.cpp -I$(INCLUDE) -o $(OBJ)/gft_superpixels.o


$(OBJ)/gft_evaluation.o: $(SRC)/gft_evaluation.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_evaluation.cpp -I$(INCLUDE) -o $(OBJ)/gft_evaluation.o

$(OBJ)/gft_marker.o: $(SRC)/gft_marker.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_marker.cpp -I$(INCLUDE) -o $(OBJ)/gft_marker.o


$(OBJ)/gft_scene.o: $(SRC)/gft_scene.cpp $(INCLUDE)/gft_scene.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_scene.cpp -I$(INCLUDE) -o $(OBJ)/gft_scene.o

$(OBJ)/gft_scene32.o: $(SRC)/gft_scene32.cpp $(INCLUDE)/gft_scene32.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_scene32.cpp -I$(INCLUDE) -o $(OBJ)/gft_scene32.o

$(OBJ)/gft_scene64.o: $(SRC)/gft_scene64.cpp $(INCLUDE)/gft_scene64.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_scene64.cpp -I$(INCLUDE) -o $(OBJ)/gft_scene64.o

$(OBJ)/gft_scene32f.o: $(SRC)/gft_scene32f.cpp $(INCLUDE)/gft_scene32f.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_scene32f.cpp -I$(INCLUDE) -o $(OBJ)/gft_scene32f.o

$(OBJ)/gft_scene64f.o: $(SRC)/gft_scene64f.cpp $(INCLUDE)/gft_scene64f.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_scene64f.cpp -I$(INCLUDE) -o $(OBJ)/gft_scene64f.o

$(OBJ)/gft_scene16.o: $(SRC)/gft_scene16.cpp $(INCLUDE)/gft_scene16.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_scene16.cpp -I$(INCLUDE) -o $(OBJ)/gft_scene16.o

$(OBJ)/gft_scene8.o: $(SRC)/gft_scene8.cpp $(INCLUDE)/gft_scene8.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_scene8.cpp -I$(INCLUDE) -o $(OBJ)/gft_scene8.o

$(OBJ)/gft_adjrel3.o: $(SRC)/gft_adjrel3.cpp $(INCLUDE)/gft_adjrel3.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_adjrel3.cpp -I$(INCLUDE) -o $(OBJ)/gft_adjrel3.o

$(OBJ)/gft_scnmath.o: $(SRC)/gft_scnmath.cpp $(INCLUDE)/gft_scnmath.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_scnmath.cpp -I$(INCLUDE) -o $(OBJ)/gft_scnmath.o

$(OBJ)/gft_matrix.o: $(SRC)/gft_matrix.cpp $(INCLUDE)/gft_matrix.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_matrix.cpp -I$(INCLUDE) -o $(OBJ)/gft_matrix.o

$(OBJ)/gft_bmap.o: $(SRC)/gft_bmap.cpp $(INCLUDE)/gft_bmap.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_bmap.cpp -I$(INCLUDE) -o $(OBJ)/gft_bmap.o

$(OBJ)/gft_set.o: $(SRC)/gft_set.cpp $(INCLUDE)/gft_set.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_set.cpp -I$(INCLUDE) -o $(OBJ)/gft_set.o

$(OBJ)/gft_stack.o: $(SRC)/gft_stack.cpp $(INCLUDE)/gft_stack.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_stack.cpp -I$(INCLUDE) -o $(OBJ)/gft_stack.o

$(OBJ)/gft_analysis.o: $(SRC)/gft_analysis.cpp $(INCLUDE)/gft_analysis.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_analysis.cpp -I$(INCLUDE) -o $(OBJ)/gft_analysis.o

$(OBJ)/gft_pqueue32.o: $(SRC)/gft_pqueue32.cpp $(INCLUDE)/gft_pqueue32.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_pqueue32.cpp -I$(INCLUDE) -o $(OBJ)/gft_pqueue32.o

$(OBJ)/gft_pqueue16.o: $(SRC)/gft_pqueue16.cpp $(INCLUDE)/gft_pqueue16.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_pqueue16.cpp -I$(INCLUDE) -o $(OBJ)/gft_pqueue16.o

$(OBJ)/gft_hashtable.o: $(SRC)/gft_hashtable.cpp $(INCLUDE)/gft_hashtable.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_hashtable.cpp -I$(INCLUDE) -o $(OBJ)/gft_hashtable.o

$(OBJ)/gft_curve.o: $(SRC)/gft_curve.cpp $(INCLUDE)/gft_curve.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_curve.cpp -I$(INCLUDE) -o $(OBJ)/gft_curve.o

$(OBJ)/gft_radiometric.o: $(SRC)/gft_radiometric.cpp $(INCLUDE)/gft_radiometric.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_radiometric.cpp -I$(INCLUDE) -o $(OBJ)/gft_radiometric.o

$(OBJ)/gft_radiometric3.o: $(SRC)/gft_radiometric3.cpp $(INCLUDE)/gft_radiometric3.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_radiometric3.cpp -I$(INCLUDE) -o $(OBJ)/gft_radiometric3.o

$(OBJ)/gft_filtering.o: $(SRC)/gft_filtering.cpp $(INCLUDE)/gft_filtering.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_filtering.cpp -I$(INCLUDE) -o $(OBJ)/gft_filtering.o

$(OBJ)/gft_filtering3.o: $(SRC)/gft_filtering3.cpp $(INCLUDE)/gft_filtering3.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_filtering3.cpp -I$(INCLUDE) -o $(OBJ)/gft_filtering3.o

$(OBJ)/gft_adjregion3.o: $(SRC)/gft_adjregion3.cpp $(INCLUDE)/gft_adjregion3.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_adjregion3.cpp -I$(INCLUDE) -o $(OBJ)/gft_adjregion3.o

$(OBJ)/gft_seedmap3.o: $(SRC)/gft_seedmap3.cpp $(INCLUDE)/gft_seedmap3.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_seedmap3.cpp -I$(INCLUDE) -o $(OBJ)/gft_seedmap3.o

$(OBJ)/gft_segmentation3.o: $(SRC)/gft_segmentation3.cpp $(INCLUDE)/gft_segmentation3.h
	$(CXX) $(FLAGS) -c $(SRC)/gft_segmentation3.cpp -I$(INCLUDE) -o $(OBJ)/gft_segmentation3.o

$(OBJ)/gft_features.o: $(SRC)/gft_features.cpp
	$(CXX) $(FLAGS) -c $(SRC)/gft_features.cpp -I$(INCLUDE) -o $(OBJ)/gft_features.o

$(OBJ)/gft_bzlib.o: $(SRC)/gft_bzlib.c
	gcc $(FLAGS) -c $(SRC)/gft_bzlib.c -I$(INCLUDE) -o $(OBJ)/gft_bzlib.o

$(OBJ)/znzlib.o: $(SRC)/znzlib.c
	gcc $(FLAGS) -c $(SRC)/znzlib.c -I$(INCLUDE) -o $(OBJ)/znzlib.o 

$(OBJ)/nifti1_io.o: $(SRC)/nifti1_io.c
	gcc $(FLAGS) -c $(SRC)/nifti1_io.c -I$(INCLUDE) -o $(OBJ)/nifti1_io.o



$(BIN)/superpixels: $(DEMO)/superpixels.cpp libgft
	$(CXX) $(FLAGS) $(DEMO)/superpixels.cpp $(GFTFLAGS) $(GFTLIB) -o $(BIN)/superpixels $(LINKS)

$(BIN)/supervoxels: $(DEMO)/supervoxels.cpp libgft
	$(CXX) $(FLAGS) $(DEMO)/supervoxels.cpp $(GFTFLAGS) $(GFTLIB) -o $(BIN)/supervoxels $(LINKS)

clean:
	rm -rf $(LIB); rm -rf $(OBJ); rm -rf $(BIN);


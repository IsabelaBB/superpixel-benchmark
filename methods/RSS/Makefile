CXX = g++ -g
CXXFLAGS = -Wall -fPIC -pedantic -Wno-unused-result -O3 -fopenmp -O3 -lstdc++ -ffast-math -march=skylake -mfma -Wno-unused-result
LIBS = -lm `pkg-config --cflags --libs opencv4`


all:
	$(CXX) $(CXXFLAGS) -c RSS.cpp -o RSS.o $(LIBS)
	$(CXX) $(CXXFLAGS) main.cpp RSSV.cpp RSS.o -o RSS_demo RSS.h $(LIBS)

clean:
	rm -f RSS.o
	rm -f RSS_demo

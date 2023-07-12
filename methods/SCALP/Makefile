CFLAGS=-Wall 
LIBS=-O2 -L/usr/X11R6/lib -lm -lpthread -lX11
EXE=SCALP

all:
	g++ $(CFLAGS) SCALP.cpp -o $(EXE) $(LIBS)

clean:
	rm -f SCALP;

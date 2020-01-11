CC = g++

ROOTLIB=$(shell root-config --libs) -lMathMore
ROOTINC=$(shell root-config --cflags)

FLAG =  -O3 -g -D_FILE_OFFSET_BITS=64
MFILE_LIB = 
RUNFILE = gsort

INCLUDE =
LIBS = 

sourcefile=main.cpp gsort.cpp gsort.h

all: $(RUNFILE)

$(RUNFILE): $(sourcefile)
	$(CC) $(FLAG) $(ROOTLIB) $(ROOTINC) -o $@ $(filter %.cpp, $(sourcefile))

clean:
	rm -f $(RUNFILE)

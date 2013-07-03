# Makefile
# CFLAGS = -m64 -O3
CFLAGS = -c -O3
LFLAGS = -lz -lgsl -lgslcblas -llapack -lboost_filesystem 
# -lgzstream -lblas -lg2c 
EXECUTABLE = altra
EXECUTABLE_2 = sim_sam
LIBS = -L/usr/local/lib/
HEADER = src/utils.h
HEADER_2 = src/Arguments.h
# -L/mnt/lustre/home/epantaleo/src/GZSTREAM/gzstream
# INCLUDES = -I/mnt/lustre/home/epantaleo/src/GZSTREAM/gzstream

${EXECUTABLE}: ${HEADER} ${HEADER_2} src/${EXECUTABLE}.cpp
	g++ -g src/${EXECUTABLE}.cpp -o $@ ${LFLAGS} ${LIBS} 

clean:
	rm -f *~ ${OBJS} ${MOREOBJS}

${EXECUTABLE_2}: ${HEADER} src/${EXECUTABLE_2}.cpp
	g++ -g src/${EXECUTABLE_2}.cpp  -o $@ ${LFLAGS} ${LIBS}

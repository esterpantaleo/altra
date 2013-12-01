SHELL = /bin/sh
LFLAGS = -lz -lgsl -lgslcblas -llapack -lboost_filesystem 
LIBS = -L/usr/local/lib/
OBJS = src/*.o
CC = g++
CFLAGS = #-Wall -Wextra
BIN = $HOME/bin
OBJs = src/altra.o src/expression.o src/TranscriptModel.o src/likelihood.o src/Reads.o src/ProposedTranscript.o src/Move.o src/Centers.o src/Transcript.o src/Center.o src/Arguments.o src/utils.o

all: ${OBJs}
	make altra

altra: ${OBJs}
	$(CC) $(CFLAGS) ${OBJs} -o src/altra ${LFLAGS} ${LIBS}

src/altra.o: src/altra.cpp src/expression.o
	$(CC) $(CFLAGS) -c $< -o $@

src/expression.o: src/expression.cpp src/TranscriptModel.o
	$(CC) $(CFLAGS) -c $< -o $@

src/TranscriptModel.o: src/TranscriptModel.cpp src/ProposedTranscript.o src/Reads.o src/Move.o src/likelihood.o 
	$(CC) $(CFLAGS) -c $< -o $@

src/likelihood.o: src/likelihood.cpp src/utils.o src/Arguments.o
	$(CC) $(CFLAGS) -c $< -o $@

src/Reads.o: src/Reads.cpp src/Transcript.o
	$(CC) $(CFLAGS) -c $< -o $@

src/ProposedTranscript.o: src/ProposedTranscript.cpp src/Centers.o 
	$(CC) $(CFLAGS) -c $< -o $@

src/Move.o: src/Move.cpp src/Centers.o
	$(CC) $(CFLAGS) -c $< -o $@

src/Centers.o: src/Centers.cpp src/Transcript.o 
	$(CC) $(CFLAGS) -c $< -o $@

src/Transcript.o: src/Transcript.cpp src/Center.o
	$(CC) $(CFLAGS) -c $< -o $@

src/Center.o: src/Center.cpp src/utils.o src/Arguments.o
	$(CC) $(CFLAGS) -c $< -o $@

src/Arguments.o: src/Arguments.cpp src/utils.o
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *~ ${OBJS} 

src/utils.o: src/utils.cpp
	$(CC) -c $< -o $@

#install: altra
#	cp altra $(BIN)

check: all
	./test/test.sh >& test/test.log
	@echo "see test/test.log"

# Makefile
SHELL = /bin/sh
LFLAGS = -lz -lgsl -lgslcblas -llapack -lboost_filesystem 
LIBS = -L/usr/local/lib/
OBJS = src/*.o
CC = g++

altra: src/altra.o
	$(CC) src/altra.o src/expression.o src/TranscriptModel.o src/likelihood.o src/Reads.o src/ProposedTranscript.o src/Move.o src/Centers.o src/Transcript.o src/Center.o src/Arguments.o src/utils.o -o $@ ${LFLAGS} ${LIBS} 
	rm -f *~ ${OBJS}

src/altra.o: src/altra.cpp src/expression.o
	$(CC) -c $< -o $@

src/expression.o: src/expression.cpp src/TranscriptModel.o
	$(CC) -c $< -o $@

src/TranscriptModel.o: src/TranscriptModel.cpp src/ProposedTranscript.o src/Reads.o src/Move.o src/likelihood.o 
	$(CC) -c $< -o $@

src/likelihood.o: src/likelihood.cpp src/utils.o src/Arguments.o
	$(CC) -c $< -o $@

src/Reads.o: src/Reads.cpp src/Transcript.o
	$(CC) -c $< -o $@

src/ProposedTranscript.o: src/ProposedTranscript.cpp src/Centers.o 
	$(CC) -c $< -o $@

src/Move.o: src/Move.cpp src/Centers.o
	$(CC) -c $< -o $@

src/Centers.o: src/Centers.cpp src/Transcript.o 
	$(CC) -c $< -o $@

src/Transcript.o: src/Transcript.cpp src/Center.o
	$(CC) -c $< -o $@

src/Center.o: src/Center.cpp src/utils.o src/Arguments.o
	$(CC) -c $< -o $@

src/Arguments.o: src/Arguments.cpp src/utils.o
	$(CC) -c $< -o $@

clean:
	rm -f *~ ${OBJS} 

sim_sam: src/sim_sam.o src/utils.o
	$(CC) $^ -o $@ ${LFLAGS} ${LIBS}

src/sim_sam: src/sim_sam.cpp src/utils.o
	$(CC) -c $< -o $@ ${LFLAGS} ${LIBS}
	rm -f *~ ${OBJS}

src/utils.o: src/utils.cpp
	$(CC) -c $< -o $@

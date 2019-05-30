# Makefile for compiling VCFtoSummStats

TARGET = VCFtoSummStats
CC = g++

# the default target that gets built when you type 'make':
all: ${TARGET}

# rule for build:
${TARGET}: ${TARGET}.cpp
	${CC} ${TARGET}.cpp -o ${TARGET}

# rule for cleaning up everything:
clean:
	rm -f ${TARGET}



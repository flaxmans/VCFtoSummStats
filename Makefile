# Makefile for compiling VCFtoSummStats

TARGET = VCFtoSummStats
CC = g++
CCFLAGS = -g

# the default target that gets built when you type 'make':
all: ${TARGET}

# rule for build:
${TARGET}: ${TARGET}.cpp ${TARGET}.hpp
	${CC} ${TARGET}.cpp ${CCFLAGS} -o ${TARGET}

# rule for cleaning up everything:
clean:
	rm -f ${TARGET}



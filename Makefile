# Makefile for compiling VCFtoSummStats

TARGET = VCFtoSummStats
CC = g++
LFLAGS = -lboost_iostreams

# conditional compiling:
DEBUG_MODE?=n
ifeq "$(DEBUG_MODE)" "Y" 
        CCFLAGS=-g -DDEBUG
else 
        CCFLAGS=-O3
endif

# the default target that gets built when you type 'make':
all: ${TARGET}

# rule for build:
${TARGET}: ${TARGET}.cpp ${TARGET}.hpp
	${CC} ${CCFLAGS} ${TARGET}.cpp ${LFLAGS} -o ${TARGET}

# rule for cleaning up everything:
clean:
	rm -f ${TARGET}



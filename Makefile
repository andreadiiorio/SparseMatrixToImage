#Copyright Andrea Di Iorio 2021
#This file is part of sparseMatrixToImage
#sparseMatrixToImage is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#sparseMatrixToImage is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with sparseMatrixToImage.  If not, see <http://www.gnu.org/licenses/>.


CC=gcc
CFLAGS=-Wall -Wextra -O2 
#libs
CFLAGS+=-lm -fopenmp
#TODO extra to reduce useless warnings
CFLAGS+=-Wno-pointer-sign -Wno-unused-parameter -Wno-unused-but-set-variable
CFLAGS+=-Wno-unused-label -Wno-unused-function
##TODO OMP CONFIG
#OMP_CANCELLATION=true
#export OMP_CANCELLATION

#SYSTEM CONFIGURATION
UNAME=$(shell uname -a | tr -c -d \[:alnum:\] | tr \[:lower:\] \[:upper:\] ) #upper uname-a
TMPDIR=/run/user/$(shell id -u)/

CONSTS = -DTMPDIR='"$(TMPDIR)"'     #TODO PUT IN /raccolta

MACROS = -DDEBUGPRINT="if(FALSE)" -DDEBUG="if(FALSE)" -DCONSISTENCY_CHECKS="if(FALSE)" -DVERBOSE="if(FALSE)"  -DDEBUGCHECKS="if(FALSE)"
MACROSDBG = -DCONSISTENCY_CHECKS="if(TRUE)"  -DDEBUGCHECKS="if(TRUE)" -DVERBOSE="if(TRUE)" -DDEBUG="if(TRUE)"
UNDEF := $(shell echo $(MACROSDBG) | tr " " "\n" | grep -oe '-D.*=' | tr -d "=" |  sed s/-D/-U/ )

CFLAGS+=$(RUNTIME)
#bind to source original project

objs := $(shell  grep -Eo '.*\..*:\s' Makefile | grep -v -e '@' -e PHONY | awk -F: '{print $1}' | tr '\n:' ' ' )
all: $(objs)

sparseMatrixToImage.o:	sparseMatrixToImage.c commons/*.c lib/*.c include/*.h 
	$(CC) $(CFLAGS) -Iinclude -I. $(filter-out %.h,$^) -o $@ $(CONSTS) $(MACROS) \
	-D MAIN_SPMAT_IMG 
sparseMatrixToImageTest.o:	sparseMatrixToImage.c commons/*.c lib/*.c include/*.h
	$(CC) $(CFLAGS) -Iinclude -I. $(filter-out %.h,$^) -o $@ $(CONSTS) $(MACROS) \
		-D MAIN_SPMAT_IMG -DTEST -ggdb -O0
clean:
	rm -f *.o

try:   
	echo $(objs) 
	echo $(UNDEF)

.PHONY: all clean testAll

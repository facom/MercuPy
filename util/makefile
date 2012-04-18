################################################################################
# MercuPy Wrap
# Jorge Zuluaga (C) 2011
################################################################################
ARCH=$(shell uname -r | cut -f 3 -d '-')
CC=g++
CFLAGS_OPT=-O4
CSPICE_PREFIX=cspice
CPPFLAGS=$(CFLAGS_OPT) -I. -I$(CSPICE_PREFIX)/include
CSPICE_LIBDIR=$(CSPICE_PREFIX)/lib
UTILPROGS=elem2state.out state2elem.out

ifeq ($(ARCH),amd64)
	CSPICE_LIBDIR=$(CSPICE_PREFIX)/lib64
	EXT=e64
endif
ifeq ($(ARCH),generic)
	CSPICE_LIBDIR=$(CSPICE_PREFIX)/lib64
	EXT=e64
endif
CSPICE_LIBS=$(CSPICE_LIBDIR)/cspice.a $(CSPICE_LIBDIR)/csupport.a
LFLAGS=-lm $(CSPICE_LIBS) 

all:$(UTILPROGS)
	@echo "Compiling..." $^

install:
	cp -rf $(UTILPROGS:.out=) ../bin

%.out:%.o
	@echo "Compiling..." $@
	${CC} $< ${LFLAGS} -o $@
	cp -rf $@ $(@:.out=)
	mv $@ $(@:out=$(EXT))

clean:
	@echo "Cleaning..."
	rm -rf *.out *.e* *~
	rm -rf $(UTILPROGS:.out=)
	cd ../bin;rm -rf $(UTILPROGS:.out=)

showarch:
	@echo $(ARCH)

################################################################################
# MercuPy Wrap
# Jorge Zuluaga (C) 2011
################################################################################
FC=gfortran
PROGRAMS=mercury6_2.exe element6.exe close6.exe
FFLAGS=-w
UTIL_PREFIX=
ifdef CSPICE
	UTIL_PREFIX=CSPICE=$(CSPICE)
endif

#################################################################################
#BASIC RULS
#################################################################################
all:$(PROGRAMS)

%.exe:%.for
	$(FC) $(FFLAGS) $^ -o $@ &> errors.log

utilbuild:
	$(UTIL_PREFIX) make -s -C util all install

updaterepo:
	@echo "Commiting..."
	@git commit -am "Commit"
	@echo "Pushing..."
	@git push origin master

pull:
	@echo "Pulling..."
	@git pull

#################################################################################
#MERCUPY PIPELINE
#################################################################################
prepare:
	bin/mercupy-prepare

build:
	make all

run:$(PROGRAMS)
	@sleep 1
	time ./mercury6_2.exe
	$(shell \
	if [ -e "output/phase" ];\
	then cp -rf output/phase output/prevphase;\
	fi)
	./element6.sh
#THIS LINE WAS ADDED FOR THE RUN PIPELINE OF SCI2WEB
	date +%s.%N > end.sig

close:
	./close6.sh

out:
	bin/mercupy-output

ref:
	bin/mercupy-ref2ref

plot:
	bin/mercupy-plot

continue:
	bin/mercupy-continue

go:
	make run out ref plot continue

err:
	bin/mercupy-diff output/$(BODY).dat.PH$(PH1) output/$(BODY).dat.PH$(PH2) 
	BODY_NAME=$(BODY) gnuplot bin/bodyerrors.gpl

errspice:
	bin/mercupy-diff SSI/$(BODY)_SPICE.dat output/$(BODY).dat 
	BODY_NAME=$(BODY)_SPICE gnuplot bin/bodyerrors.gpl

#################################################################################
#CLEAN RULES
#################################################################################
clean:cleanprepare

cleanprepare:cleanbuild
	rm -rf  *.in *.inc

cleanbuild:cleanrun
	rm -rf *.exe

cleanrun:cleanout
	rm -rf *.out *.dmp *.tmp
	rm -rf *.state
	rm -rf *.oxt
	rm -rf fort.*

cleanout:cleanref
	rm -rf output/*.*
	rm -rf output/*phase*

cleanref:cleanplot
	rm -rf output/*.ref

cleanplot:cleantmp
	rm -rf output/*.png

cleantmp:
	rm -rf *~ *.tmp *.dmp *.log tmp/*
	find . -name "*~" -exec rm -rf {} \;

#ADDED FOR THE RUN PIPELINE OF SCI2WEB
cleansci2web:
	rm -rf *.sig
	rm -rf run.sh run.conf submit.sh

commit:
	git commit -am "Commit" && git push origin master

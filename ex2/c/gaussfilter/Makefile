# Compiler
CC = gcc
CFLAGS += -Wall -Werror -Wextra -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

# For gnuplot
SCRIPT = log.gpi
FORMAT = eps
PLOTTER = gnuplot
void =
TAB = $(\t)

TARGETS = myPlot.eps output.txt input.txt
DATAFILE = output.txt


DATAIMPORT = ./../../Diffraction/Si100_180sec.mca
START = $(shell sed -n '/<<DATA>>/ =' $(DATAIMPORT))
END   = $(shell sed -n '/<<END>>/  =' $(DATAIMPORT))
#==============================================================================#


NUMBEROFPEAKS = 1

#LINESTART = $(shell for i in `seq 1 1 $(NUMBEROFPEAKS)`; do echo `sed -n "/<<START$$i>>/ =" $(DATAFILE)`)

LINEMARKERS = $(shell sed -n '/<<START0>>/ =' $(DATAFILE)) 
LINESTART1  = $(shell sed -n '/<<START1>>/ =' $(DATAFILE))
LINESTART2  = $(shell sed -n '/<<START2>>/ =' $(DATAFILE))
LINESTART = $(LINESTART0) $(LINESTART1) $(LINESTART2)

LINEEND0    = $(shell sed -n '/<<END0>>/ =' $(DATAFILE))
LINEEND1    = $(shell sed -n '/<<END1>>/ =' $(DATAFILE))
LINEEND2    = $(shell sed -n '/<<END2>>/ =' $(DATAFILE))
LINEEND = $(LINEEND0) $(LINEEND1) $(LINEEND2)

GAUSS0START = $(shell awk 'NR==$(LINESTART0)+1' $(DATAFILE))
GAUSS1START = $(shell awk 'NR==$(LINESTART1)+1' $(DATAFILE))
GAUSS2START = $(shell awk 'NR==$(LINESTART2)+1' $(DATAFILE))
GAUSSSTART = $(GAUSS0START) $(GAUSS1START) $(GAUSS2START)


GAUSS0END   = $(shell awk 'NR==$(LINEEND0)-1' $(DATAFILE))
GAUSS1END   = $(shell awk 'NR==$(LINEEND1)-1' $(DATAFILE))
GAUSS2END   = $(shell awk 'NR==$(LINEEND2)-1' $(DATAFILE))
GAUSSEND = $(GAUSS0END) $(GAUSS1END) $(GAUSS2END)

# THis written in shorter way?


#.PHONEY: all
all: $(TARGETS)

#==============================================================================#

# Figure 1
myPlot.$(FORMAT): input.txt output.txt Makefile
	echo 'set terminal eps' > $(SCRIPT)
	echo 'set output "$@"' >> $(SCRIPT)
	echo 'set key out' >> $(SCRIPT)
	echo 'set title "This is a title"' >> $(SCRIPT)
	echo 'set grid' >> $(SCRIPT)
	echo 'set xlabel "x"' >> $(SCRIPT)
	echo 'set ylabel "y"' >> $(SCRIPT)
	echo "set datafile sep ','">> $(SCRIPT)
	echo 'plot \$(void)' >> $(SCRIPT)
	echo '"input.txt" using 1:2 every ::$(GAUSS0START)::$(GAUSS0END) with points title "gauss0"\$(void)' >>$(SCRIPT)
	echo ', "input.txt" using 1:2 every ::$(GAUSS1START)::$(GAUSS1END) with points title "gauss1"\$(void)' >>$(SCRIPT)
	echo ', "input.txt" using 1:2 every ::$(GAUSS2START)::$(GAUSS2END) with points title "gauss2"' >>$(SCRIPT)


	$(PLOTTER) $(SCRIPT)
	gv $@
#==============================================================================#

# Output of main
output.txt: main input.txt
	./$< $(NUMBEROFPEAKS) $(word 2, $^) > $@

#==============================================================================#

# Input
input.txt: $(DATALFILE)
	$(shell awk 'NR==$(START)+1, NR==$(END)-1' $^ > input.txt)

.PHONEY:clean
clean:                     # this is "clean" target. it does not depend on anything
	find ./ -type f -executable -delete
	$(RM) *.dat *.svg *.png *.pdf log* *.eps



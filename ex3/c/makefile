CC = gcc
CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

# For gnuplot
SCRIPT = log.gpi
FORMAT = eps
PLOTTER = gnuplot

# Remember:
# $^ (prerequisite)
# $@ (target)

#.PHONEY: all
all: main

# Output of main
output.txt: main
	./$^ > $@

# Implicit rules to produce main from .c file
main: main.c


.PHONEY: clean
clean :
	find ./ -type f -executable -delete
	$(RM) *.dat $(SCRIPT) *.svg *.png *.pdf log* *.eps

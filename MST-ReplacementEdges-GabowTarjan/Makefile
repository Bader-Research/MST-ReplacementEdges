CC = gcc
PROG = mst
FILES := MST-GabowTarjan.c disjointsets-GabowTarjan.c
GRAPHS = $(wildcard inputs/*.gr)
OUTFILES = $(patsubst inputs/%.gr, %.out, $(GRAPHS))

.PHONY: test clean tarball

$(PROG) : $(FILES)
#	$(CC) -O3 -o $(PROG) $(FILES) -lrt
	$(CC) -O3 -g -Wall -o $(PROG) $(FILES)

%.out : inputs/%.gr $(PROG)
	-./$(PROG) -q -f $< -o $@

test: $(PROG) $(OUTFILES)

tarball:
	tar cvf $(PROG).tar $(FILES) $(patsubst %.c,%.h,$(FILES)) $(GRAPHS) Makefile
	gzip $(PROG).tar

clean :
	-rm -f *.o $(PROG)

DCP = dcp
DCP_SOLVER = dcp_solver
BSGS = bsgs
BSGS_SOLVER = bsgs_solver
DCP_GLV_SOLVER = dcp_glv_solver
TEST_PARI = test_pari
SHELL  = /bin/sh
CFLAGS     = -O3 -Wall -ffp-contract=off -fno-strict-aliasing
CC         = /usr/bin/gcc
CPPFLAGS   = -I. -I/usr/local/include
LD         = /usr/bin/gcc
LDFLAGS    = -O3 -Wall -ffp-contract=off -fno-strict-aliasing    -Wl,--export-dynamic 
RUNPTH     = -Wl,-rpath "/usr/local/lib"
DLCFLAGS   = -fPIC
LIBS       = -lm -L/usr/local/lib -lpari
RM = rm -f
ALL = $(DCP) $(DCP_SOLVER) $(BSGS) $(DCP_GLV_SOLVER) $(BSGS_SOLVER)


all: $(ALL)

$(DCP): $(DCP).c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<

$(DCP_SOLVER): $(DCP_SOLVER).c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
	$(LD) -o $@ $(LDFLAGS) $(DCP_SOLVER).o $(DCP).o $(RUNPTH) $(LIBS)
	
$(BSGS): $(BSGS).cpp
	/usr/bin/g++ -c $(CFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
	mkdir -p results

$(DCP_GLV_SOLVER): $(DCP_GLV_SOLVER).c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
	$(LD) -o $@ $(LDFLAGS) $(DCP_GLV_SOLVER).o $(DCP).o $(RUNPTH) $(LIBS)

$(BSGS_SOLVER): $(BSGS_SOLVER).cpp
	/usr/bin/g++ -c $(CFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
	/usr/bin/g++ -o $@ $(LDFLAGS) $(BSGS_SOLVER).o $(BSGS).o  $(RUNPTH) $(LIBS)


test: $(TEST_PARI).cpp $(BSGS).cpp
	/usr/bin/g++ -c $(CFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
	/usr/bin/g++ -o $@ $(LDFLAGS) $(TEST_PARI).o $(BSGS).o $(RUNPTH) $(LIBS)	

clean:
	-$(RM) *.o $(ALL) test
	-$(RM) -r results/*

# Irons (1970) frontal solution program: CODEST, ZIPP, ELMAK + tests
FC       ?= gfortran
F77FLAGS  = -std=legacy -O2
F90FLAGS  = -O2
# ZIPP's error exit deliberately divides by zero (Irons' trace artifice);
# -ffpe-summary=none stops gfortran printing an IEEE note about it.
LDFLAGS   = -ffpe-summary=none

all: test_irons_deck test_irons_fe

test_irons_deck: zcapture.o zipp.o elmak_card.o test_irons_deck.o
	$(FC) $(LDFLAGS) -o $@ $^
test_irons_fe: zcapture.o zipp.o test_irons_fe.o
	$(FC) $(LDFLAGS) -o $@ $^

test_irons_deck.o test_irons_fe.o: zcapture.o

%.o: %.f
	$(FC) $(F77FLAGS) -c $<
%.o: %.f90
	$(FC) $(F90FLAGS) -c $<

check: all
	./test_irons_deck | tee irons_results.txt
	./test_irons_fe  | tee -a irons_results.txt

clean:
	rm -f *.o *.mod test_irons_deck test_irons_fe irons_*.dat irons_*.lis irons_results.txt

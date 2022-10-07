FC = gfortran
FCFLAGS = -fbacktrace -fno-align-commons -fbounds-check -std=legacy
FLFLAGS = -lblas



SRCS = main1.f Hwang1d.f vode.f dgbfa.f dgbsl.f dgefa.f dgesl.f  

OBJS = main1.o Hwang1d.o vode.o dgbfa.o dgbsl.o dgefa.o dgesl.o 
#OBJS = $(patsubst %.f, %.o, $(SRCS))

PROGRAM = example

all: $(OBJS)
	$(FC) $(FCFLAGS) -o $(PROGRAM) $^ $(FLFLAGS)


$(OBJS): $(SRCS)
	$(FC) $(FCFLAGS) -c -o $@ $< $(FLFLAGS)

%.o : %.f
	$(FC) $(FCFLAGS) -c $< $(FLFLAGS)

clean: 
	rm -f *.o *.mod

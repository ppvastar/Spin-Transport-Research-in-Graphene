# $Id: Makefile,v 1.5 2002/04/0Q4 03:42:49 weng Exp $

SOURCE0 = constant.f90
SOURCE =        scat.f90   generaterou.f90	 drift.f90	initial.f90		
SOURCE += 	rkf.f90		main.f90         generateEl.f90			
SOURCE += 	coherent.f90 	generategama.f90
SOURCE += 	generateg_v_u.f90	

exec=gtransport.out

OBJ0 = $(patsubst %.f90,%.o,$(SOURCE0))
OBJ = $(patsubst %.f90,%.o,$(SOURCE))
OBJ1 = $(patsubst %.f90,%.o,$(SOURCE1))

F77 = /usr/local/intel/fce/10.1.022/bin/ifort -static -openmp
F90 = /usr/local/intel/fce/10.1.022/bin/ifort -static -openmp
FFLAGS = -O -g 
FFLAGS_FPP =  
FFLAGS_FPP += -fpp 

FFLAGS_FPP += -I/usr/local/mpich-intel/include/f90base
LDFLAGS += -L/usr/local/mpich-intel/lib
FLIBS += -lmpichf90 -lmpich

#--------------------------------------------------------------
#--------------------------------------------------------------

#LDFLAGS += -L/usr/local/intel/mkl/10.0.5.025/lib/em64t
#FLIBS += -lmkl_lapack -lmkl_em64t -lmkl_solver -lguide -lpthread
#FFLAGS := -static -openmp


#.SURFFIXES:

%.o : %.F90
	$(F90) $(FFLAGS) $(MACRO) $(FFLAGS_FPP) -c -o $@ $<
%.o : %.f90
	$(F90) $(FFLAGS) $(MACRO) $(FFLAGS_FPP) -c -o $@ $<

all:	$(exec)

$(OBJ) : $(OBJ0)

gtransport.out:	$(OBJ) $(OBJ1) $(OBJ0)
	$(F90) $(LDFLAGS) -o $@ $^ $(FLIBS)


etags:TAGS
TAGS:	$(SOURCE) $(SOURCE1) $(SOURCE0)
	etags $^

clean:	
	rm -f *.od *.o  *.mod *.dat
	rm -f $(exec) nohup.out
  

distclean:
	rm -f *.od *.o *.mod *~ 
	rm -f $(exec)
	rm -f *.dat nohup.out PI*

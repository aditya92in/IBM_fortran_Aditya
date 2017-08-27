CF=gfortran
LD=gfortran
#CF=f90
#LD=f90
#FFLAGS1=-c -Oscalar3,vector3,task3,aggress -r2  -ef -l perf
#FFLAGS2=-c -Oscalar3,vector3,task3,aggress -r2  -ef -l perf
#LDFLAGS= -ef -l perf
#FFLAGS1=-c -Oscalar3,vector3,task3,aggress -r2   -N132
#FFLAGS2=-c -Oscalar3,vector3,task3,aggress -r2   -N132
#FFLAGS1=-c -s -r8
#FFLAGS2=-c -s -r8
#FFLAGS1=-c -Rabc -N132
#FFLAGS2=-c -Oscalar3,vector3,task3,aggress -r2   -N132
#FFLAGS1=-c  -O +autodblpad +es -C 
#FFLAGS2=-c  -O +autodblpad  +es
# Checks: -Hasu
#
# -fpe0: abort on all exceptions EXCEPT underflow
#   checking: -C  -O1 -traceback

#FFLAGS1=-c  -W0 -cm  -w90  -C  -O1 -traceback        -r8 -extend_source      
FFLAGS1=-c -ffixed-form -ffixed-line-length-none #-W0 -cm -O3 -w90          -r8 -extend_source      
FFLAGS2=-c -ffixed-form -ffixed-line-length-none #-cm -O3 -w90  -r8 -extend_source 
FFLAGSNAG=-c -ffixed-form -ffixed-line-length-none #-cm -O3 -w90  -r8 
#FFLAGSBOUND=-c -W0 -cm -O3 -w90  -r8 -extend_source 
LDFLAGS=  
FFLAGS1=-c -ffixed-form -ffixed-line-length-none #-W0 -cm -O3 -w90          -r8 -extend_source      
FFLAGS2=-c -ffixed-form -ffixed-line-length-none #-cm -O3 -w90  -r8 -extend_source 
FFLAGSNAG=-c -ffixed-form -ffixed-line-length-none #-cm -O3 -w90  -r8 
#FFLAGS1=-c  -fw -w -Ad -Kfast -O3    -Wa,-W     
#FFLAGS2=-c  -Ad  -fw  -w -Ad -Kfast -O3
#FFLAGSNAG=-c  -w  -Kfast -O3 -X9  
#FFLAGSbound=-c  -Ad -w  -Kfast -O3   
#LDFLAGS= -static 
#LIB = -lnag_intel
EXE=model.e
OBJSNAG=g05cbf.o g05daf.o timer.o
OBJS1=\
MAIN.o fillps.o loadd.o      \
fkdat.o  chkdiv.o tec.o     reagr2.o   \
chkdtm.o mea2.o    adams.o   poisson.o  \
corgen.o boucha.o bouchp.o momve5.o\
adamsles.o momles.o \
setp0.o \
bouio.o boupio.o \
tec_cropped.o  \
adscal.o bouioc.o\
stepscal.o 



#FFT = fft.o 
FFT = fft.o cycred.o
LIBDIR=lib
#LIB=$(LIBDIR)/$(EXE).a

$(EXE):		$(OBJS1) $(OBJSf) $(FFT)  $(OBJSNAG)
		$(LD)  -o $@ $(OBJS1) $(OBJSf) $(FFT)  $(OBJSNAG) $(LDFLAGS) $(LIB) $(SLIBS) 

$(OBJS1):	%.o: %.f param.txt common.txt  
		$(CF) $(FFLAGS1) $(SFLAGS) $<


$(FFT):		%.o: %.f
		$(CF) $(FFLAGS2) $(SFLAGS) $<

$(OBJSNAG):	%.o: %.f 
		$(CF) $(FFLAGSNAG) $(SFLAGS) $<


clean:
		rm -rf $(EXE) *.[ol] core
		touch core
		chmod -w core

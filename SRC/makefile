

pre_flags   = ifort -O2 -xHost

post_flags  = -L${LD_LIBRARY_PATH} -I${MKLROOT}/include -qopenmp -lpthread
final_flags = -L${LD_LIBRARY_PATH} -I${MKLROOT}/include -Wl,--start-group  ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread

my_sources=$(wildcard *.f)
my_objects=$(my_sources:.f=.o)

other_objects = libxdrf.a 

treecode3d_final.o: treecode3d_final.f
	$(pre_flags) -c treecode3d_final.f $(post_flags)

uiowa_bd.modules.o:      uiowa_bd.modules.f
	$(pre_flags) -c uiowa_bd.modules.f $(post_flags) 

uiowa_bd.o:      uiowa_bd.f
	$(pre_flags) -c uiowa_bd.f $(post_flags) 

# the following should deal with all subroutine files here:
%.o: %.f 
	$(pre_flags) -c $*.f -o $*.o $(post_flags)

uiowa_bd: $(my_objects) $(other_objects)
	$(pre_flags) $(my_objects) $(other_objects) -o uiowa_bd.exe $(final_flags)

clean:
	rm *.mod
	rm *.o

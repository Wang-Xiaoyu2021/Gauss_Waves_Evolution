gauss : gauss.f90
	gfortran -o gauss gauss.f90 disturbance.f90
	./gauss
.PHONY : clean
clean : 
	rm gauss realoutput.dat
	make clean



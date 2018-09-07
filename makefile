gauss : gauss.f90
	gfortran -o gauss gauss.f90 inv_mat.f90
	./gauss
.PHONY : clean
clean : 
	rm gauss realoutput.dat
	make clean



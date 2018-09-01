gauss : gauss.f90
	gfortran -o gauss gauss.f90
	./gauss
.PHONY : clean
clean : 
	rm gauss
	make clean



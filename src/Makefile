FLAGS=
FC=gfortran

FC_W_OPTS=$(FC) $(FLAGS)
LINKER=$(FC) $(FLAGS)

../gradient_descent_exp: lin_reg_fort.o rotate_integrate.o gradient_descent_exp.o 
	$(LINKER) rotate_integrate.o lin_reg_fort.o gradient_descent_exp.o -o ../gradient_descent_exp

lin_reg_fort.o: lin_reg_fort.f90
	$(FC_W_OPTS) -c lin_reg_fort.f90
rotate_integrate.o: rotate_integrate.f90
	$(FC_W_OPTS) -c rotate_integrate.f90
gradient_descent_exp.o: gradient_descent_exp.f90
	$(FC_W_OPTS) -c gradient_descent_exp.f90

clean:
	rm -f *.o *.mod ../gradient_descent_exp

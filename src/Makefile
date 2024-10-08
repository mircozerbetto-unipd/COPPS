#######################################################################################################################################
##                                                  CONFIGURATION PART                                                               ##
#######################################################################################################################################

include ../Makefile.in

#######################################################################################################################################
##                                                     BUILDING PART                                                                 ##
#######################################################################################################################################
copps: basis.o bessel.o copps.o experimental.o euler.o lanczos.o matrix.o observables.o physics.o powell.o relax.o rq.o s3j.o spline.o stvec.o tensor.o tensor4.o tensor5.o wigner.o
	$(CC) $(FLAG) $(II) $(LL) basis.o bessel.o copps.o experimental.o euler.o lanczos.o matrix.o observables.o physics.o powell.o relax.o rq.o s3j.o spline.o stvec.o tensor.o tensor4.o tensor5.cpp wigner.o -lcqp++ -lminpack_jac -lmv -lsparse -lspblas -llevmar $(LAPACKLIB) $(FORTLIB) -lcubature -ljemdi -lcmrg -lm -o copps
basis.o : basis.cpp prep.h 
	$(CC) $(II) $(LL) -c $(FLAG) basis.cpp
bessel.o : bessel.cpp
	$(CC) $(II) $(LL) -c $(FLAG) bessel.cpp
copps.o : copps.cpp copps.h prep.h
	$(CC) $(II) $(LL) -c $(FLAG) copps.cpp  
experimental.o : experimental.cpp experimental.h prep.h
	$(CC) $(II) $(LL) -c $(FLAG) experimental.cpp 
euler.o : euler.cpp euler.h
	$(CC) $(II) $(LL) -c $(FLAG) euler.cpp 
lanczos.o : lanczos.cpp lanczos.h prep.h
	$(CC) $(II) $(LL) -c $(FLAG) lanczos.cpp 
matrix.o : matrix.cpp matrix.h prep.h
	$(CC) $(II) $(LL) -c $(FLAG) matrix.cpp 
observables.o : observables.cpp copps.h
	$(CC) $(II) $(LL) -c $(FLAG) observables.cpp 
physics.o : physics.cpp physics.h prep.h
	$(CC) $(II) $(LL) -c $(FLAG) physics.cpp
powell.o : powell.f
	$(FF) $(II) $(LL) -c $(FLAGF) powell.f
relax.o : relax.cpp relax.h prep.h
	$(CC) $(II) $(LL) -c $(FLAG) relax.cpp 
rq.o : rq.cpp
	$(CC) $(II) $(LL) -c $(FLAG) rq.cpp
s3j.o : s3j.cpp s3j.h
	$(CC) $(II) $(LL) -c $(FLAG) s3j.cpp 
spline.o : spline.cpp spline.h
	$(CC) $(II) $(LL) -c $(FLAG) spline.cpp 
stvec.o : stvec.cpp stvec.h prep.h
	$(CC) $(II) $(LL) -c $(FLAG) stvec.cpp 
tensor.o : tensor.cpp tensor.h
	$(CC) $(II) $(LL) -c $(FLAG) tensor.cpp 
tensor4.o : tensor4.cpp tensor4.h
	$(CC) $(II) $(LL) -c $(FLAG) tensor4.cpp 
tensor5.o : tensor5.cpp tensor5.h
	$(CC) $(II) $(LL) -c $(FLAG) tensor5.cpp 
wigner.o : wigner.cpp wigner.h
	$(CC) $(II) $(LL) -c $(FLAG) wigner.cpp 

clean:
	rm -f ./*o ./copps




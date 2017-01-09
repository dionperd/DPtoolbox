%mex -lgsl -lgslcblas /usr/lib/libgfortran.so.3 dlsode.o integration_loop.c -o integration_loop.mexa64

%mex -lgsl -lgslcblas /usr/lib/libgfortran.so.3 dlsode.o integration_loop_I.c -o integration_loop_I.mexa64

%mex -lgsl -lgslcblas /usr/lib/libgfortran.so.3 dlsode.o integration_loop_I_EXC.c -o integration_loop_I_EXC.mexa64

%mex -lgsl -lgslcblas integration_loop.c -o integration_loop.mexa64

% mex -lgsl -lgslcblas integration_loop_seq.c -o integration_loop_seq.mexa64

%mex -lgsl -lgslcblas integration_loop_big_seq.c -o integration_loop_big_seq.mexa64

%mex -lgsl -lgslcblas integration_loop_handwriting.c -o integration_loop_handwriting.mexa64

%mex -lgsl -lgslcblas integration_loop_handwriting_lamdaNoNauto.c -o integration_loop_handwriting_lamdaNoNauto.mexa64

%mex -lgsl -lgslcblas integration_loop_handwriting_XiNoNauto.c -o integration_loop_handwriting_XiNoNauto.mexa64

%mex -lgsl -lgslcblas integration_loop_handwriting_XiNoNautoEuler.c -o integration_loop_handwriting_XiNoNautoEuler.mexa64



%for my macbook pro
mex -I/Users/denis/Documents/Software/gsl_universal_1.14/include -L/Users/denis/Documents/Software/gsl_universal_1.14/libs -lgsl -lgslcblas integration_loop_handwriting.c -o integration_loop_handwriting.mexmaci

mex -I/opt/local/include/ -L/opt/local/lib -lgsl -lgslcblas integration_loop_handwriting.c -o integration_loop_handwriting.mexmaci

mex -I/opt/local/var/macports/software/gsl/1.14_0/opt/local/include/ -L/opt/local/var/macports/software/gsl/1.14_0/opt/local/lib -lgsl -lgslcblas integration_loop_handwriting.c -o integration_loop_handwriting.mexmaci

%for my Mac-Os-X
mex -I/usr/local/include/ -L/usr/local/lib -lgsl -lgslcblas integration_loop.c -o integration_loop.mexmaci64


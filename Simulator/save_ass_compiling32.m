%mex -lgsl -lgslcblas /usr/lib/libgfortran.so.3 dlsode.o integration_loop.c -o integration_loop

%mex -lgsl -lgslcblas /usr/lib/libgfortran.so.3 dlsode.o integration_loop_I.c -o integration_loop_I

%mex -lgsl -lgslcblas /usr/lib/libgfortran.so.3 dlsode.o integration_loop_I_EXC.c -o integration_loop_I_EXC

%mex -lgsl -lgslcblas integration_loop.c -o integration_loop

% mex -lgsl -lgslcblas integration_loop_seq.c -o integration_loop_seq

%mex -lgsl -lgslcblas integration_loop_big_seq.c -o integration_loop_big_seq

mex -lgsl -lgslcblas integration_loop_handwriting.c -o integration_loop_handwriting32

% compile
disp('Making mexGmmTrainSP...');
mex mexGmmTrainSP.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
disp('Making mexGmmTrainDP...');
mex mexGmmTrainDP.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

disp('Making mexFisherEncodeHelperSP...');
mex mexFisherEncodeHelperSP.cxx ../fisher.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"

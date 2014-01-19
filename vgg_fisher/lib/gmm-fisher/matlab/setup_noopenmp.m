% compile
disp('Making mexGmmTrainSP...');
mex mexGmmTrainSP.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims
disp('Making mexGmmTrainDP...');
mex mexGmmTrainDP.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims

disp('Making mexFisherEncodeHelperSP...');
mex mexFisherEncodeHelperSP.cxx ../fisher.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims

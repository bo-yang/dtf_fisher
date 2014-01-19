% compile
disp('Making mexGmmTrainSP...');
mex -g mexGmmTrainSP.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims
disp('Making mexGmmTrainDP...');
mex -g mexGmmTrainDP.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims

disp('Making mexFisherEncodeHelperSP...');
mex -g mexFisherEncodeHelperSP.cxx ../fisher.cxx ../gmm.cxx ../stat.cxx ../simd_math.cxx -largeArrayDims

MATLAB="/opt/well/matlab/matlab-2013b"
Arch=glnxa64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/home/cyau/.matlab/R2013b"
OPTSFILE_NAME="./mexopts.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for logsumexp2" > logsumexp2_mex.mki
echo "CC=$CC" >> logsumexp2_mex.mki
echo "CFLAGS=$CFLAGS" >> logsumexp2_mex.mki
echo "CLIBS=$CLIBS" >> logsumexp2_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> logsumexp2_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> logsumexp2_mex.mki
echo "CXX=$CXX" >> logsumexp2_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> logsumexp2_mex.mki
echo "CXXLIBS=$CXXLIBS" >> logsumexp2_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> logsumexp2_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> logsumexp2_mex.mki
echo "LD=$LD" >> logsumexp2_mex.mki
echo "LDFLAGS=$LDFLAGS" >> logsumexp2_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> logsumexp2_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> logsumexp2_mex.mki
echo "Arch=$Arch" >> logsumexp2_mex.mki
echo OMPFLAGS= >> logsumexp2_mex.mki
echo OMPLINKFLAGS= >> logsumexp2_mex.mki
echo "EMC_COMPILER=" >> logsumexp2_mex.mki
echo "EMC_CONFIG=optim" >> logsumexp2_mex.mki
"/opt/well/matlab/matlab-2013b/bin/glnxa64/gmake" -B -f logsumexp2_mex.mk

TARGET="$1"
BUILD_MODE="$2"
if [[ -z "$TARGET" ]]; then
    echo "Usage: ./build.sh [cpu|gpu]"
    exit 1
fi

# PATHS
HOME_PATH="$(dirname "$PWD")"
SPACK_PATH=$HOME_PATH/spack
. $SPACK_PATH/share/spack/setup-env.sh 
PROJECT_DIR=$HOME_PATH/x3d2
BUILD_DIR="$PROJECT_DIR/build_$TARGET"

# Set flags depending on the mode
if [[ "$TARGET" == "gpu" ]]; then
    	echo "Building with GPU flags..."
	# CMAKE_FLAGS="-DUSE_GPU=ON -DUSE_CPU=OFF"
	spack env activate -p $PROJECT_DIR/scripts/spack-env-gpu
	spack install
	spack load nvhpc
	spack load cmake	
	# Fatal Error: Cannot open module file 'mpi.mod' for reading at (1): No such file or directory
	# spack location -i nvhpc	
elif [[ "$TARGET" == "cpu" ]]; then
    	echo "Building with CPU flags..."
	spack env activate -p $PROJECT_DIR/scripts/spack-env-cpu
	spack install
	spack load openmpi
	spack load cmake
else
    echo "Invalid option: $TARGET"
    echo "Use: cpu or gpu"
    exit 1
fi

export FC=mpif90
rm -rf $BUILD_DIR
if [[ "$BUILD_MODE" == "debug" ]]; then
	cmake -S "$PROJECT_DIR" -B $BUILD_DIR -DCMAKE_BUILD_TYPE=Debug
else
	cmake -S "$PROJECT_DIR" -B $BUILD_DIR -DCMAKE_BUILD_TYPE=Release
fi
cd $BUILD_DIR
make
echo "Running Tests..."
make test

set(CMAKE_SYSTEM_NAME Linux)

set(CMAKE_C_COMPILER   x86_64-conda-linux-gnu-gcc)
set(CMAKE_CXX_COMPILER x86_64-conda-linux-gnu-g++)

# where is the target environment located
set(CMAKE_FIND_ROOT_PATH  $ENV{HOME}/a1/apps/miniforge3/envs/sortmerna-build)

# find_program on host system
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)

# search headers and libraries in the target environment only
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

set(CMAKE_SYSROOT $ENV{HOME}/a1/apps/miniforge3/envs/sortmerna-build/x86_64-conda-linux-gnu/sysroot)
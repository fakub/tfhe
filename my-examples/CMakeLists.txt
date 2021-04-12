cmake_minimum_required(VERSION 3.0)

#   set 'optim' build type by default
set(CMAKE_CONFIGURATION_TYPES Debug Optim CACHE TYPE INTERNAL FORCE)
set(CMAKE_BUILD_TYPE "optim" CACHE STRING "Build Type: Debug or Optim")

#   project name
project(my-examples-tfhe)

#   C/C++ flags
set(CMAKE_CXX_FLAGS_DEBUG "${CLANG_FLAGS} -std=gnu++11 -g3 -O0 -Wall -Werror")
set(CMAKE_C_FLAGS_DEBUG "-g3 -O0 -Wall -Werror")

set(CMAKE_CXX_FLAGS_OPTIM "${CLANG_FLAGS} -std=gnu++11 -g3 -march=native -O2 -DNDEBUG -funroll-loops -Wall -Werror")
set(CMAKE_C_FLAGS_OPTIM "-g3 -march=native -O3 -DNDEBUG -funroll-loops -Wall -Werror")

set(CMAKE_CXX_FLAGS_RELEASE "${CLANG_FLAGS} -std=gnu++11 -g0 -march=native -O2 -DNDEBUG -funroll-loops -Wall -Werror")
set(CMAKE_C_FLAGS_RELEASE "-g0 -march=native -O3 -DNDEBUG -funroll-loops -Wall -Werror")

#   include & lib dirs
include_directories("../bin/include/tfhe")
file(GLOB TFHE_HEADERS ../bin/include/tfhe/*.h)
link_directories("../bin/lib")
set(RUNTIME_LIBS tfhe-nayuki-portable)

#   list of examples to build
set(CPP_EXAMPLES
        example-orig-addition-boot
        example-parallel-addition
        )

#   build an executable for each example
foreach (CPP_EXAMPLE ${CPP_EXAMPLES})
    add_executable(${CPP_EXAMPLE}-nayuki-portable ${CPP_EXAMPLE}.cpp ${TFHE_HEADERS})
    target_link_libraries(${CPP_EXAMPLE}-nayuki-portable ${RUNTIME_LIBS})
endforeach (CPP_EXAMPLE)

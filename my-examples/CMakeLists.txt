cmake_minimum_required(VERSION 3.0)

project(my-examples-tfhe)

#~ # -std=c99 seems to be required in Travis tests for whatever reason
#~ set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c99")

set(CMAKE_CXX_FLAGS_DEBUG "${CLANG_FLAGS} -std=gnu++11 -g3 -O0 -Wall -Werror")
set(CMAKE_C_FLAGS_DEBUG "-g3 -O0 -Wall -Werror")

set(CMAKE_CXX_FLAGS_OPTIM "${CLANG_FLAGS} -std=gnu++11 -g3 -march=native -O2 -DNDEBUG -funroll-loops -Wall -Werror")
set(CMAKE_C_FLAGS_OPTIM "-g3 -march=native -O3 -DNDEBUG -funroll-loops -Wall -Werror")

set(CMAKE_CXX_FLAGS_RELEASE "${CLANG_FLAGS} -std=gnu++11 -g0 -march=native -O2 -DNDEBUG -funroll-loops -Wall -Werror")
set(CMAKE_C_FLAGS_RELEASE "-g0 -march=native -O3 -DNDEBUG -funroll-loops -Wall -Werror")

include_directories("../bin/include/tfhe")
file(GLOB TFHE_HEADERS ../bin/include/tfhe/*.h)
link_directories("../bin/lib")


set(CPP_EXAMPLES
        example-addition-boot
        )


set(RUNTIME_LIBS tfhe-nayuki-portable)

#the integration tests must be single source code, and are compiled as a standalone application
#we first compile the C++ tests
foreach (CPP_EXAMPLE ${CPP_EXAMPLES})
    add_executable(${CPP_EXAMPLE}-nayuki-portable ${CPP_EXAMPLE}.cpp ${TFHE_HEADERS})
    target_link_libraries(${CPP_EXAMPLE}-nayuki-portable ${RUNTIME_LIBS})
endforeach (CPP_EXAMPLE)

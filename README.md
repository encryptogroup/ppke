# Efficient and Privacy-Preserving Kidney Exchange Protocol <br> SPIKE: Secure and Private Investigation of the Kidney Exchange problem
This is part of a joint work of Timm Birka, Kay Hamacher, Tobias Kussel, Helen Möllering, and Thomas Schneider.

This work is based on the [ABY framework](https://github.com/encryptogroup/ABY/) for efficient mixed-protocol secure two-party computation.

For simplicity, we will only describe the short version of the build process here, including the differences and additional steps required compared to the ABY build procedure. Please also refer to the detailed build instructions of ABY.

This code is provided as an experimental implementation for testing purposes and should not be used in a productive environment. We cannot guarantee security and correctness.

# Build Requirements

This code requires all [ABY requirements](https://github.com/encryptogroup/ABY#requirements) and can be set up and executed like any other [ABY example](https://github.com/encryptogroup/ABY#aby-applications). Additionally, we use the C++ [Json Library of Niels Lohmann](https://github.com/nlohmann/json) to process Json files. The Library can be included as an external dependency in the ABY CMake File by adding the following lines to the CMakeLists.txt:

```
find_package(nlohmann_json QUIET)
if(nlohmann_json_FOUND)
    message(STATUS "Found nlohmann_json")
elseif(NOT nlohmann_json_FOUND AND NOT TARGET nlohmann_json::nlohmann_json)
    message("nlohmann_json was not found: add nlohmann_json subdirectory")
    if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/nlohmann_json/CMakeLists.txt")
        find_package(Git REQUIRED)
        message("initialize Git submodule: extern/nlohmann_json")
        execute_process(COMMAND git submodule update --init extern/nlohmann_json
                        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
    endif()
    add_subdirectory(extern/nlohmann_json)
endif()
```

# Parameters Compatibility Matching
```
 -r [Role: 0/1, required]
 -a [IP-address, default: localhost, optional]
 -d [Path pointing to the location where the data is stored. Default: '../data/input/data_1000', optional]
 -f [Factors to be included 0/1, default: 1, optional]
 -x [Number of participating pairs, default: 10, optional]
 -y [Number of HLA used, default 50, optional]
```
More (optional) parameters can be inspected by simply running ./matching.
# Parameters Cycle Computation
```
 -r [Role: 0/1, required]
 -a [IP-address, default: localhost, optional]
 -x [Number of participating pairs, default: 10, optional]
 -y [Length of cycles, default: 3, optional]
```
More (optional) parameters can be inspected by simply running ./compute_cycles.
# Parameters Cycle and Solution Evaluation
```
 -r [Role: 0/1, required]
 -a [IP-address, default: localhost, optional]
 -x [Number of participating pairs, default: 10, optional]
 -y [Length of cycles, default: 3, optional]
```
More (optional) parameters can be inspected by simply running ./find_largest_set.
# Parameters Full Protocol
```
 -r [Role: 0/1, required]
 -a [IP-address, default: localhost, optional]
 -d [Path pointing to the location where the data is stored. Default: '../data/input/data_1000', optional]
 -f [Factors to be included 0/1, default: 1, optional]
 -x [Number of participating pairs. Default: 10, optional]
 -y [Number of HLA used. Default 50, optional]
 -z [Length of cycles. Default: 3, optional]
```
More (optional) parameters can be inspected by simply running ./eppkep.

# Running
We have included simple scripts ([Compatibility Matching](https://github.com/encryptogroup/ppke/tree/main/scripts/matching), [Cycle Computation](https://github.com/encryptogroup/ppke/tree/main/scripts/compute_cycles), [Cycle and Solution Evaluation](https://github.com/encryptogroup/ppke/tree/main/scripts/find_largest_set), and [EPPKEP](https://github.com/encryptogroup/ppke/tree/main/scripts/eppkep)) to test our implementation locally. Our Protocol consists of four parts (the third and fouth part are both in one [file]()). The parts can be executed independently but each part depends on the input of the previous part(s), e.g., the second part (Cycle Computation) requires the compatibility graph computed in part one, matching, the third and fourth part (Cycle and Solution Evaluation). Example data is given [here](https://github.com/encryptogroup/ppke/tree/main/data_1000). The data was generated using [python scripts](https://github.com/encryptogroup/ppke/tree/main/src/generate_test_data).

# Issues & Notes
ABY has a known problem where it randomly returns invalid results in about 1 out of 10 executions. This problem is still an open and also reported in several issues, e.g., [here](https://github.com/MPC-SoK/frameworks/issues/19) or [here](https://github.com/encryptogroup/ABY/issues/114).

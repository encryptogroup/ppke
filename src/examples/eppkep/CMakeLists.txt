
add_executable(matching matching.cpp common/matching_circuit.cpp common/patient.cpp common/pair.cpp common/plain_data.cpp)
target_link_libraries(matching ABY::aby ENCRYPTO_utils::encrypto_utils)


add_executable(compute_cycles compute_cycles.cpp common/compute_cycles_circuit.cpp)
target_link_libraries(compute_cycles ABY::aby ENCRYPTO_utils::encrypto_utils)

add_executable(find_largest_set find_largest_set.cpp common/find_largest_set_circuit.cpp)
target_link_libraries(find_largest_set ABY::aby ENCRYPTO_utils::encrypto_utils)

add_executable(eppkep eppkep.cpp 
                    common/matching_circuit.cpp common/patient.cpp common/pair.cpp common/plain_data.cpp
                    common/compute_cycles_circuit.cpp
                    common/find_largest_set_circuit.cpp)
target_link_libraries(eppkep ABY::aby ENCRYPTO_utils::encrypto_utils)
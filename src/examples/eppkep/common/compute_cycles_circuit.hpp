#ifndef __COMPUTE_CYCLES_CIRCUIT__HPP_
#define __COMPUTE_CYCLES_CIRCUIT__HPP_

#include "CircuitWrapper.h"
#include "../../../abycore/circuit/booleancircuits.h"
#include "../../../abycore/circuit/arithmeticcircuits.h"
#include "../../../abycore/circuit/circuit.h"
#include "../../../abycore/aby/abyparty.h"
#include <math.h>
#include <cassert>
#include <memory>
#include <vector>

using cmp_graph = std::vector<std::vector<share_p>>;


/**
 \param		role 		    role played by the program which can be server or client part.
 \param 	address 	    IP Address
 \param 	seclvl 		    Security level
 \param 	bitlen		    Bit length of the inputs
 \param 	nthreads	    Number of threads
 \param		mt_alg		    The algorithm for generation of multiplication triples
 \param         n_pairs             Number of pairs participating
 \param         cycle_length        Length of the cycles   
 \brief		This function is used for running a testing environment for the second part of the protocol
 */
uint32_t compute_cycles_circuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length);


/**
 * @brief Builds the circuit to compute the number of cycles within a compatibility graph.
 * 
 * @param role role played by the program which can be server or client part.
 * @param address IP Address
 * @param port Port
 * @param seclvl Security level 
 * @param bitlen Bit length of the inputs
 * @param nthreads Number of threads 
 * @param mt_alg The algorithm for generation of multiplication triples 
 * @param n_pairs Number of pairs participating
 * @param cycle_length Length of the cycles 
 * @return Number of cycles within the graph
 */
uint32_t BuildComputeCyclesCircuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length);


/**
 * @brief Reads the shared compatibility graph from files and returns it.
 * 
 * @param ac Arithmetic circuitwrapper
 * @param role  Either SERVER or CLIENT
 * @param bitlen Length of the values
 * @param n_pairs Number of participating pairs
 * @return std::vector<std::vector<share_p>> 
 */
cmp_graph readCompGraphFromFile(CircuitW_p ac, e_role role, uint32_t bitlen, uint32_t n_pairs);


/**
 * @brief Constructs an unweighted compatibility graph out of the given, weighted compatibility graph.
 * 
 * @param s_comp_graph Compatibility graph in boolean sharing
 * @param ac Arithmetic circuitwrapper
 * @param bc Boolean circuitwrapper
 * @param bitlen Length of the input values
 * @param n_pairs Number of pairs
 * @return Unweighted compatibility graph 
 */
cmp_graph BuildComputeUnweightedGraphCircuit(cmp_graph s_comp_graph, CircuitW_p ac, 
                CircuitW_p bc, uint32_t bitlen, uint32_t n_pairs);
 

/**
 * @brief Computes the cycle_length-power of the unweighted compatibility graph.
 * 
 * @param s_unweighted_graph Unweighted compatibility graph in arithmetic sharing
 * @param ac Arithmetic circuitwrapper
 * @param bc Boolean circuitwrapper
 * @param yc Yao circuitwrapper
 * @param bitlen Lenght of the input values 
 * @param n_pairs Number of pairs
 * @param cycle_length Length of cycles
 * @return Cycle_length power of the unweighted compatibility graph. 
 */
cmp_graph BuildMatrixMultiplicationCircuit(cmp_graph s_unweighted_graph, CircuitW_p ac, 
                CircuitW_p bc, CircuitW_p yc, uint32_t bitlen, uint32_t n_pairs, 
                uint32_t cycle_length);


/**
 * @brief Counts the number of cycles (entries on the diagonal) in the compatibility graph for a given length.
 * 
 * @param s_length_matrix Cycle_length-power of the unweighted compatibility graph in arithmetic sharing
 * @param ac Arithmetic circuitwrapper
 * @param bitlen Length of the input values
 * @param n_pairs Number of pairs
 * @return Sum of the values on the diagonal.
 */
share_p BuildCountCyclesCircuit(cmp_graph s_length_matrix, CircuitW_p ac, uint32_t bitlen, uint32_t n_pairs);



// Plain implementation
std::vector<std::vector<uint32_t>> readCompGraphFromFilePlain(uint32_t n_pairs);


uint32_t compute_number_of_cycles(uint32_t n_pairs, uint32_t cycle_length);

#endif 
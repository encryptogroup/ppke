#ifndef __FIND_LARGEST_SET_CIRCUIT__HPP_
#define __FIND_LARGEST_SET_CIRCUIT__HPP_

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
using cycles = std::vector<std::tuple<share_p, std::vector<share_p>>>;
using cycles_set = std::vector<std::vector<share_p>>;
using result_set = std::tuple<share_p, cycles_set>; 

/**
 \param		role 		    role played by the program which can be server or client part.
 \param 	address 	    IP Address
 \param 	seclvl 		    Security level
 \param 	bitlen		    Bit length of the inputs
 \param 	nthreads	    Number of threads
 \param		mt_alg		    The algorithm for generation of multiplication triples
 \param     n_pairs         Number of pairs participating
 \param     cycle_length    Length of the cycles   
 \brief		This function is used for running a testing environment for the second part of the protocol
 */
void find_largest_set_circuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, 
                uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length, uint32_t n_cycles);


/**
 * @brief 
 * 
 * @param role 
 * @param address 
 * @param port 
 * @param seclvl 
 * @param bitlen 
 * @param nthreads 
 * @param mt_alg 
 * @param n_pairs 
 * @param cycle_length 
 * @param n_all_cycles 
 * @param n_cycles 
 */
void findCyclesCircuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, 
                uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length, uint32_t n_all_cycles,
                uint32_t n_cycles);


/**
 * @brief 
 * 
 * @param role 
 * @param address 
 * @param port 
 * @param seclvl 
 * @param bitlen 
 * @param nthreads 
 * @param mt_alg 
 * @param n_pairs 
 * @param cycle_length 
 * @param n_all_cycles 
 * @param n_cycles 
 */
void findSetCircuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, 
                uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length, uint32_t n_cycles);

                
/**
 * @brief Reads the shared compatibility graph from files and returns it
 * 
 * @param ac Arithmetic circuitwrapper
 * @param role Either SERVER or CLIENT
 * @param bitlen Length of the values
 * @param n_pairs Number of participating pairs
 * @return cmp_graph 
 */
cmp_graph readCompGraphFromFile1(CircuitW_p &ac, e_role role, uint32_t bitlen, uint32_t n_pairs);


/**
 * @brief 
 * 
 * @param cycleA 
 * @param cycleB 
 * @param yc 
 * @param bitlen 
 * @param cycle_length 
 * @return share_p 
 */
share_p compare_cycles(std::vector<share_p> &cycleA, std::vector<share_p> &cycleB, CircuitW_p &yc, 
                uint32_t bitlen, uint32_t cycle_length);


/**
 * @brief 
 * 
 * @param all_cycles 
 * @param cycle 
 * @param weight
 * @param yc 
 * @param cycle_length 
 * @return share_p 
 */
share_p containsCycle(cycles &all_cycles, std::vector<share_p> &cycle, share_p &weight, CircuitW_p &yc, 
                uint32_t bitlen, uint32_t cycle_length);


/**
 * @brief Modifies all_cycles in place and adds every cycle it finds to it
 * 
 * @param all_cycles 
 * @param s_cmp_graph 
 * @param s_cmp_graph_gmw 
 * @param current_cycle 
 * @param current_length 
 * @param weight 
 * @param valid 
 * @param yc 
 * @param bitlen 
 * @param n_pairs 
 * @param cycle_length 
 */
void findCycles(cycles &all_cycles, cmp_graph &s_cmp_graph, std::vector<uint32_t> &current_cycle, 
                uint32_t current_length, share_p &weight, share_p &valid, CircuitW_p &yc, 
                uint32_t bitlen, uint32_t n_pairs, uint32_t cycle_length);


/**
 * @brief Removes all duplicates from the set of cycles
 * 
 * @param set_cycles 
 * @param yc 
 * @param bitlen 
 * @param n_cycles 
 * @param cycle_length 
 * @return cycles 
 */
cycles removeDuplicates(cycles &sorted_cycles, CircuitW_p &yc, uint32_t bitlen, uint32_t n_pairs, uint32_t n_cycles, 
                    uint32_t cycle_length);


/**
 * @brief 
 * 
 * @param all_cycles 
 * @param yc 
 * @param bitlen 
 * @param n_pairs 
 * @param n_all_cycles 
 * @param n_cycles 
 * @param cycle_length 
 * @return cycles 
 */
cycles kNNSort(cycles &all_cycles, CircuitW_p &yc, uint32_t bitlen, uint32_t n_pairs, uint32_t n_all_cycles, 
                uint32_t n_cycles, uint32_t cycle_length);


/**
 * @brief 
 * 
 * @param sorted_cycles 
 * @param bc 
 * @param n_cycles 
 * @param cycle_length 
 */
void prepareWritingCyclesToFile(cycles &sorted_cycles, CircuitW_p &bc, uint32_t n_cycles, uint32_t cycle_length);


/**
 * @brief 
 * 
 * @param sorted_cycles 
 * @param role 
 * @param n_cycles 
 * @param cycle_length 
 */
void writeCyclesToFile(cycles &sorted_cycles, e_role role, uint32_t n_cycles, uint32_t cycle_length);


/**
 * @brief 
 * 
 * @param bc 
 * @param role 
 * @return cycles 
 */
cycles readCyclesFromFile(CircuitW_p &bc, e_role role, uint32_t bitlen, uint32_t n_cycles, uint32_t cycle_length);


/**
 * @brief 
 * 
 * @param set_cycles 
 * @param cycle 
 * @param bc 
 * @param bitlen 
 * @param cycle_count 
 * @param cycle_length 
 * @return share_p 
 */
share_p disjointSet(cycles_set &set_cycles, std::vector<share_p> &cycle, CircuitW_p &bc, uint32_t bitlen, 
                uint32_t cycle_count, uint32_t cycle_length);


/**
 * @brief returns the set with the highest weight, implemented using a slight variation of the knn sort
 * 
 * @param all_cycle_sets 
 * @param set_weights 
 * @param yc 
 * @param bitlen 
 * @param n_cycles 
 * @param cycle_length 
 * @return result 
 */
result_set findMaximumSet(std::vector<cycles_set> &all_cycle_sets, std::vector<share_p> &set_weights, 
                CircuitW_p &yc, uint32_t bitlen, uint32_t n_pairs, uint32_t n_cycles, 
                uint32_t cycle_length);


/**
 * @brief 
 * 
 * @param cycles 
 * @param bc 
 * @param yc 
 * @param bitlen 
 * @param n_cycles 
 * @param cycle_length 
 * @return std::vector<std::vector<share_p>> 
 */
result_set findSolution(cycles &all_cycles, CircuitW_p &bc, CircuitW_p &yc, uint32_t bitlen, 
                uint32_t n_pairs, uint32_t n_cycles, uint32_t cycle_length);



// Plain implementation
void plain_solution(uint32_t n_pairs, uint32_t cycle_length, uint32_t n_all_cycles, uint32_t n_cycles);


std::vector<std::vector<uint32_t>> readCompGraphFromFilePlain1(uint32_t n_pairs);


std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> sortCyclesPlain(std::vector<std::tuple<uint32_t, 
                std::vector<uint32_t>>> cycles, uint32_t n_all_cycles);


bool containsCyclePlain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> &all_cycles, std::vector<uint32_t> &cycle, 
            uint32_t weight, uint32_t cycle_length);


std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> findCyclesPlain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> all_cycles,
            std::vector<std::vector<uint32_t>> cmp_graph, std::vector<uint32_t> current_cycle, uint32_t current_length, uint32_t weight, 
            uint32_t valid, uint32_t n_pairs, uint32_t cycle_length);


bool disjointSetPlain(std::vector<std::vector<uint32_t>> cycle_set, std::vector<uint32_t> cycle, 
                uint32_t cycle_count, uint32_t cycle_length);


std::tuple<uint32_t, std::vector<std::vector<uint32_t>>> findMaximumSetPlain(std::vector<std::vector<std::vector<uint32_t>>> cycle_sets, 
                std::vector<uint32_t> set_weights, uint32_t n_cycles);


std::tuple<uint32_t, std::vector<std::vector<uint32_t>>> findSolutionPlain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> cycles, 
                uint32_t n_pairs, uint32_t n_cycles, uint32_t cycle_length);
#endif 
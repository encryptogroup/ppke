#ifndef __MATCHING_CIRCUIT__HPP_
#define __MATCHING_CIRCUIT__HPP_

#include "CircuitWrapper.h"
#include "../../../abycore/circuit/booleancircuits.h"
#include "../../../abycore/circuit/arithmeticcircuits.h"
#include "../../../abycore/circuit/circuit.h"
#include "../../../abycore/aby/abyparty.h"
#include <math.h>
#include <cassert>
#include <memory>
#include <tuple>
#include <vector>

class Pair;
class PlainPair;

using cmp_graph = std::vector<std::vector<share_p>>;

/**
 \param		role 		role played by the program which can be server or client part.
 \param 	address 	IP Address
 \param 	seclvl 		Security level
 \param 	bitlen		Bit length of the inputs
 \param 	nthreads	Number of threads
 \param		mt_alg		The algorithm for generation of multiplication triples
 \param         path_to_data    Path pointing to the directory where the test data is stored
 \param         n_pairs         Number of pairs participating in the protocol
 \brief		This function is used for running a testing environment for the first part of the protocol
 */
void matching_circuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, std::string path_to_data, 
                const bool factors, uint32_t n_pairs, uint32_t n_hla);


/**
 * @brief 
 * 
 * @param role Either SERVER or CLIENT 
 * @param address IP adress
 * @param port Port
 * @param seclvl Security level
 * @param bitlen Length of the input values
 * @param nthreads Number of threads
 * @param mt_alg The algorithm for generation of multiplication triples
 * @param path_to_data Path pointing to the directory where the test data is stored
 * @param factors Sets whether or not all factors should be used
 * @param n_pairs Number of pairs
 * @param n_hla Number of observed hla
 */
void BuildMatchingCircuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, std::string path_to_data, 
                const bool factors, uint32_t n_pairs, uint32_t n_hla);


/**
 * @brief Reads the pair data from file and stores it in secret sharings.
 * 
 * @param bc Boolean Circuitwrapper
 * @param role Either Server or Client
 * @param bitlen Length of the inputs
 * @param path Path to data
 * @param n_pairs Number of pairs
 * @return Pair data in boolean sharing. 
 */
std::vector<std::shared_ptr<Pair>> read_pair_data(CircuitW_p bc, e_role role, uint32_t bitlen, 
                std::string &path, uint32_t n_pairs);


/**
 * @brief Performs HLA crossmatch.
 * 
 * @param s_hla_d  HLA of the donor as SIMD in boolean sharing
 * @param s_hla_a_r  HLA antibodies of the recipient as SIMD in boolean sharing
 * @param bc Boolean circuitwrapper
 * @param n_hla Number of HLA
 * @return Result of the HLA crossmatch.
 */
share_p BuildCrossmatchCircuit(share_p s_hla_d, share_p s_hla_a_r, CircuitW_p bc, uint32_t n_hla);


/**
 * @brief Performs HLA match
 * 
 * @param s_hla_d HLA of the donor as SIMD in boolean sharing
 * @param s_hla_r HLA of the recipient as SIMD in boolean sharing
 * @param s_best Constant for the best result
 * @param s_good Constant for the second best result
 * @param s_ok Constant for the third best result
 * @param s_zero Constant for the worst result
 * @param bc BooleanCircuitWrapper
 * @return Result of the HLA match
 */
share_p BuildHLAmatchCircuit(share_p s_hla_d, share_p s_hla_r, share_p s_best, share_p s_good, 
                share_p s_ok, share_p s_zero, CircuitW_p bc);


/**
 * @brief Performs ABO match.
 * 
 * @param s_abo_d ABO of the donor in boolean sharing
 * @param s_abo_r ABO of the recipient in boolean sharing
 * @param s_best Constant for the best result
 * @param s_zero Constant for the worst result
 * @param bc Boolean Circuitwrapper
 * @return Result of the ABO match.
 */
share_p BuildABOmatchCircuit(std::vector<share_p> &s_abo_d, std::vector<share_p> &s_abo_r, 
                share_p s_best, share_p s_zero, CircuitW_p bc);


/**
 * @brief Performs Age match.
 * 
 * @param s_age_d Age group of the donor in boolean sharing
 * @param s_age_r Age group of the recipient in boolean sharing
 * @param s_best Constant for the best result
 * @param s_good Constant for the second best result
 * @param s_zero Constant for the worst result
 * @param bc Boolean circuitwrapper
 * @return Result of the age match.
 */
share_p BuildAgematchCircuit(share_p s_age_d, share_p s_age_r, share_p s_best, share_p s_good, 
                share_p s_zero, CircuitW_p bc);


/**
 * @brief Performs Sex match.
 * 
 * @param s_sex_d Sex of the donor as boolean sharing
 * @param s_sex_r Sex of the recipient as boolean sharing
 * @param s_best Constant for the best result
 * @param s_good Constant for the second best result
 * @param s_zero Constant for the worst result
 * @param bc Boolean circuitwrapper
 * @return Result of the sex match.
 */
share_p BuildSexmatchCircuit(share_p s_sex_d, share_p s_sex_r, share_p s_best, share_p s_good, 
                share_p s_zero, CircuitW_p bc);


/**
 * @brief Performs size match.
 * 
 * @param s_size_d Size of the donor as boolean sharing
 * @param s_size_r Size of the recipient as boolean sharing
 * @param s_good Constant for the best result
 * @param s_zero Constant for the worst result
 * @param bc Boolean circuitwrapper
 * @return Result of the size match.
 */
share_p BuildSizematchCircuit(share_p s_size_d, share_p s_size_r, share_p s_good, share_p s_zero, 
                CircuitW_p bc);


/**
 * @brief Prepares compatibility graph for writing it to a file
 * 
 * @param comp_graph Compatibility graph in arithmetic sharing
 * @param ac Arithmetic circuitwrapper
 * @param n_pairs Number of pairs
 */
void prepareWritingCompGraphToFile(cmp_graph &comp_graph, CircuitW_p ac,
                uint32_t n_pairs);


/**
 * @brief Writes the compatibility to file while being secretshared.
 * 
 * @param comp_graph Compatibility graph in arithmetic sharing
 * @param role Either SERVER or CLIENT
 * @param n_pairs Number of pairs
 */
void writeCompGraphToFile(cmp_graph &comp_graph, e_role role, uint32_t n_pairs);



//Plain implementation
void testing(std::string path_to_data, const bool factors, uint32_t n_pairs, uint32_t n_hla);


std::vector<std::shared_ptr<PlainPair>> read_pair_data_plain(std::string &path, uint32_t n_pairs); 


uint32_t evalCrossmatch(uint32_t *hla_d, uint32_t *hla_a_r, uint32_t n_hla);


uint32_t evalHLAmatch(uint32_t *hla_d, uint32_t *hla_r, uint32_t n_hla);


uint32_t evalABOmatch(uint32_t *abo_d, uint32_t *abo_r);


uint32_t evalAgematch(uint32_t age_d, uint32_t age_r);


uint32_t evalSexmatch(uint32_t sex_d, uint32_t sex_r);


uint32_t evalSizematch(uint32_t size_d, uint32_t size_r);

void writeCompGraphToFilePlain(std::vector<std::vector<uint32_t>> &comp_graph, uint32_t n_pairs);

#endif 

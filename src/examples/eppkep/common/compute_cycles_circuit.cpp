#include "compute_cycles_circuit.hpp"
#include "../../../abycore/circuit/booleancircuits.h"
#include "../../../abycore/sharing/sharing.h"

// Json Lib
#include "../../../extern/nlohmann_json/single_include/nlohmann/json.hpp"

#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using json = nlohmann::json;

uint32_t compute_cycles_circuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		        uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length) 
{
    auto result = BuildComputeCyclesCircuit(role, address, port, seclvl, bitlen, nthreads, mt_alg, n_pairs, cycle_length);
    
    auto result_plain = compute_number_of_cycles(n_pairs, cycle_length);

    std::cout << "\n\n#####\t#####\t#####\t#####\t#####" << std::endl;
    std::cout << "Number of cycles:\t" << (int) result  << "\tExpected:\t" << (int) result_plain<< std::endl;
    std::cout << "Number of unique cycles:\t" << (int) result / (int) cycle_length << "\tExpected:\t" << (int) result_plain / (int) cycle_length << std::endl;
    // std::cout << "Number of cycles without duplicates:\t" << (int) result / (int) cycle_length << "\tExpected:\t" << (int) result_plain / (int) cycle_length << std::endl;
    // std::cout << "Maximumnumber of unique cycles (no shared vertices):\t" << (int) n_pairs / (int) cycle_length << std::endl;
    
    return result;   
}


uint32_t BuildComputeCyclesCircuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		        uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length) 
{
    ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
			mt_alg);

	std::vector<Sharing*>& sharings = party->GetSharings();

    // Create circuits -> Main computation is done with arithmetic shares
	CircuitW_p arithmcirc = std::make_shared<CircuitWrapper>(sharings[S_ARITH]->GetCircuitBuildRoutine());
    CircuitW_p yaocirc = std::make_shared<CircuitWrapper>(sharings[S_YAO]->GetCircuitBuildRoutine());
    CircuitW_p boolcirc = std::make_shared<CircuitWrapper>(sharings[S_BOOL]->GetCircuitBuildRoutine());

    auto s_comp_graph = readCompGraphFromFile(arithmcirc, role, bitlen, n_pairs); 

    // std::cout << "Read compatibility graph from file.\n" << std::endl;

    // Convert sharing to boolean sharing A2B
    for(size_t i = 0; i < n_pairs; ++i) {
        for(size_t j = 0; j < n_pairs; ++j) {
            s_comp_graph[i][j] = boolcirc->PutA2BGate(s_comp_graph[i][j], yaocirc->circ_);
        }
    }

    // std::cout << "Converted sharing from arithmetic to boolean.\n" << std::endl;

    auto s_unweighted_graph = BuildComputeUnweightedGraphCircuit(s_comp_graph, arithmcirc, boolcirc, bitlen, n_pairs);

    // std::cout << "Computed unweighted adjecency matrix.\n" << std::endl;

    auto s_length_matrix = BuildMatrixMultiplicationCircuit(s_unweighted_graph, arithmcirc, boolcirc, yaocirc, bitlen, n_pairs, cycle_length);

    // std::cout << "Computed the " << (int) cycle_length << " power of the unweighted adjacency matrix.\n" << std::endl;

    auto s_sum = BuildCountCyclesCircuit(s_length_matrix, arithmcirc, bitlen, n_pairs);
    
    // std::cout << "Counted the number of existing cycles.\n" << std::endl;

    s_sum = arithmcirc->PutOUTGate(s_sum, ALL);

    party->ExecCircuit();

    auto n_cycles = s_sum->get_clear_value<uint32_t>();

    // free memory
    party->Reset();
    delete party;

    return n_cycles;
}


cmp_graph readCompGraphFromFile(CircuitW_p ac, e_role role, uint32_t bitlen, uint32_t n_pairs) 
{
    cmp_graph s_cmp_graph;

    std::string roleString = "Client";
    if(role == SERVER) {
        roleString = "Server";
    }
    std::ifstream infile("../data/output/comp_graph_"+roleString+".json");

    json jomp_graph;
    infile >> jomp_graph;

    if(((uint32_t) jomp_graph["n_pairs"]) != n_pairs) {
        std::cout << "Number of pairs do not match!" <<std::endl;
        return s_cmp_graph;
    }
    
    json graph = jomp_graph["graph"];
    for(size_t i = 0; i < n_pairs; ++i) {
        std::vector<share_p> row;
        for(size_t j = 0; j < n_pairs; ++j) {
            row.push_back(ac->PutSharedINGate((uint32_t) graph[i][j], bitlen));
        }
        s_cmp_graph.push_back(row);
    }

    return s_cmp_graph;
}


cmp_graph BuildComputeUnweightedGraphCircuit(cmp_graph s_comp_graph, CircuitW_p ac, 
                CircuitW_p bc, uint32_t bitlen, uint32_t n_pairs) 
{    
    std::vector<std::vector<share_p>> s_unweighted_graph;
    share_p s_one = bc->PutCONSGate((uint32_t) 1, bitlen);
    share_p s_zero = bc->PutCONSGate((uint32_t) 0, bitlen);
    for(size_t i = 0; i < n_pairs; ++i) {
        std::vector<share_p> current_row;
        for(size_t j = 0; j < n_pairs; ++j) {
            share_p sel =  bc->PutGTGate(s_comp_graph[i][j], s_zero);
            share_p res =  bc->PutMUXGate(s_one, s_zero, sel);
            // Put the share in Arithmic sharing into the resulting matrix
            current_row.push_back(ac->PutB2AGate(res));
        }
        s_unweighted_graph.push_back(current_row);
    }

    return s_unweighted_graph;
}


cmp_graph BuildMatrixMultiplicationCircuit(cmp_graph s_unweighted_graph, CircuitW_p ac, 
                CircuitW_p bc, CircuitW_p yc, uint32_t bitlen, uint32_t n_pairs, 
                uint32_t cycle_length) 
{
    std::vector<std::vector<share_p>> s_length_matrix = s_unweighted_graph;

    // Start matrix multiplication - cycle lengths time
    for(size_t c = 1; c < cycle_length; ++c) {
        for(size_t i = 0; i < n_pairs; ++i) {
            std::vector<share_p> current_row;
            for(size_t j = 0; j < n_pairs; ++j) {
                current_row.push_back(ac->PutINGate((uint32_t) 0, bitlen, SERVER));
            }
            for(size_t j = 0; j < n_pairs; ++j) {
                for(size_t k = 0; k < n_pairs; ++k) {
                    share_p tmp_mul = ac->PutMULGate(s_length_matrix[i][k], s_unweighted_graph[k][j]);
                    current_row[j] = ac->PutADDGate(current_row[j], tmp_mul);
                }
            }
            s_length_matrix[i] = current_row;
        }    
    }

    return s_length_matrix;
}


share_p BuildCountCyclesCircuit(cmp_graph s_length_matrix, CircuitW_p ac, uint32_t bitlen,
                uint32_t n_pairs) 
{    
    share_p s_sum = ac->PutINGate((uint32_t) 0 , bitlen, SERVER);
    for(size_t i = 0; i < n_pairs; ++i) {
        s_sum = ac->PutADDGate(s_sum, s_length_matrix[i][i]);
    }

    return s_sum;
}



// Plain implementation
std::vector<std::vector<uint32_t>> readCompGraphFromFilePlain(uint32_t n_pairs) 
{
    std::vector<std::vector<uint32_t>> comp_graph;

    std::ifstream infile("../data/output/comp_graph_plain.json");
    std::string tmp;

    json jomp_graph;
    infile >> jomp_graph;

    if(((uint32_t)jomp_graph["n_pairs"]) != n_pairs) {
        std::cout << "Number of pairs do not match!" <<std::endl;
        return comp_graph;
    }
    
    json graph = jomp_graph["graph"];
    for(size_t i = 0; i < n_pairs; ++i) {
        std::vector<uint32_t> row;
        for(size_t j = 0; j < n_pairs; ++j) {
            row.push_back((uint32_t) graph[i][j]);
        }
        comp_graph.push_back(row);
    }

    return comp_graph;
}


uint32_t compute_number_of_cycles(uint32_t n_pairs, uint32_t cycle_length) 
{
    std::vector<std::vector<uint32_t>> comp_graph = readCompGraphFromFilePlain(n_pairs);
    
    std::vector<std::vector<uint32_t>> unweighted_graph;
    for(size_t i = 0; i < n_pairs; ++i) {
        std::vector<uint32_t> current;
        for(size_t j = 0; j < n_pairs; ++j) {
            current.push_back(comp_graph[i][j] > 0 ? 1 : 0);
        }
        unweighted_graph.push_back(current);
    }

    std::vector<std::vector<uint32_t>> resulting_matrix = unweighted_graph;
    for(size_t c = 1; c < cycle_length; ++c) {
        for(size_t i = 0; i < n_pairs; ++i) {
            std::vector<uint32_t> current(comp_graph.size());
            for(size_t j = 0; j < n_pairs; ++j) {
                for(size_t k = 0; k < n_pairs; ++k) {
                    current[j] += resulting_matrix[i][k] * unweighted_graph[k][j];
                }
            }
            resulting_matrix[i] = current;
        }
    }

    // std::cout << "#####\t#####\t#####\t#####\t#####" << std::endl;
    // for(size_t i = 0; i < n_pairs; ++i) {
    //     for(size_t j = 0; j < n_pairs; ++j) {
    //         std::cout << resulting_matrix[i][j] << " ";
    //     }
    //     std::cout << "\n" << std::endl;
    // }

    uint32_t n_cycles = 0;
    for(size_t i = 0; i < n_pairs; ++i) {
            n_cycles += resulting_matrix[i][i];
    }

    return n_cycles;
}
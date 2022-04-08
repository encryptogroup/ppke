#include "find_largest_set_circuit.hpp"
#include "../../../abycore/circuit/booleancircuits.h"
#include "../../../abycore/sharing/sharing.h"

// Json Lib
#include "../../../extern/nlohmann_json/single_include/nlohmann/json.hpp"

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using json = nlohmann::json;

uint32_t all_cycles(uint32_t n_pairs, uint32_t cycle_length) 
{
    uint32_t result = n_pairs;
    for(size_t i = 1; i < cycle_length; ++i) {
        result *= (n_pairs - i);
    }
    return result;
}


void find_largest_set_circuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		        uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, 
                uint32_t cycle_length, uint32_t n_cycles) 
{
    uint32_t n_all_cycles = all_cycles(n_pairs, cycle_length);

    findCyclesCircuit(role, address, port, seclvl, bitlen, nthreads, mt_alg, n_pairs, cycle_length, n_all_cycles,
                n_cycles);
    
    std::cout <<"\n\n##########\nFound Cycles. Continue with search for the largest set...\n##########\n" << std::endl;

    uint32_t unique_cycles = (int) n_cycles / (int) cycle_length;
    findSetCircuit(role, address, port, seclvl, bitlen, nthreads, mt_alg, n_pairs, cycle_length, unique_cycles);
    
    // std::cout << "\n#####\t#####\t#####\t#####\t#####" << std::endl;
    // std::cout << "#####\t#####\t#####\t#####\t#####\n" << std::endl;

    plain_solution(n_pairs, cycle_length, n_all_cycles, unique_cycles); 
}   


void findCyclesCircuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, 
                uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length, uint32_t n_all_cycles,
                uint32_t n_cycles)
{
    ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                mt_alg);

	std::vector<Sharing*>& sharings = party->GetSharings();

	CircuitW_p arithmcirc = std::make_shared<CircuitWrapper>(sharings[S_ARITH]->GetCircuitBuildRoutine());
    CircuitW_p yaocirc = std::make_shared<CircuitWrapper>(sharings[S_YAO]->GetCircuitBuildRoutine());
    CircuitW_p boolcirc = std::make_shared<CircuitWrapper>(sharings[S_BOOL]->GetCircuitBuildRoutine());

    auto s_cmp_graph = readCompGraphFromFile1(arithmcirc, role, bitlen, n_pairs);
    
    // Convert sharing to yao sharing
    for(size_t i = 0; i < n_pairs; ++i) {
        for(size_t j = 0; j < n_pairs; ++j) {
            s_cmp_graph[i][j] = yaocirc->PutA2YGate(s_cmp_graph[i][j]);
        }
    }

    // Prepare for find cycles
    cycles all_cycles;
    share_p weight = yaocirc->PutINGate((uint32_t) 0, bitlen, SERVER);
    share_p valid = yaocirc->PutINGate((uint32_t) 0, bitlen, SERVER);
    std::vector<uint32_t> current_cycle;
    findCycles(all_cycles, s_cmp_graph, current_cycle, 0, weight, valid, yaocirc, 
                bitlen, n_pairs, cycle_length);

    auto sorted_cycles = kNNSort(all_cycles, yaocirc, bitlen, n_pairs, n_all_cycles, n_cycles, cycle_length);

    // Remove all duplicates
    auto filtered_cycles = removeDuplicates(sorted_cycles, yaocirc, bitlen, n_pairs, n_cycles, cycle_length);
    
    // Number of unique cycles
    uint32_t unique_cycles = (int) n_cycles / (int) cycle_length;

    // Convert to GMW sharing for saving to file
    for(size_t i = 0; i < unique_cycles; ++i) {
        std::vector<share_p> current_row;
        auto cycle = std::get<1>(filtered_cycles[i]);
        for(size_t j = 0; j < cycle_length; ++j) {
            current_row.push_back(boolcirc->PutY2BGate(cycle[j]));
        }
        auto weight = boolcirc->PutY2BGate(std::get<0>(filtered_cycles[i]));
        filtered_cycles[i] = std::make_tuple(weight, current_row);
    }

    prepareWritingCyclesToFile(filtered_cycles, boolcirc, unique_cycles, cycle_length);

    party->ExecCircuit();

    writeCyclesToFile(filtered_cycles, role, unique_cycles, cycle_length);

    // Reset circuit to free memory
    party->Reset();
    delete party;    
}


void findSetCircuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, 
                uint32_t nthreads, e_mt_gen_alg mt_alg, uint32_t n_pairs, uint32_t cycle_length, uint32_t n_cycles)
{
    ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
                mt_alg);

	std::vector<Sharing*>& sharings = party->GetSharings();

    CircuitW_p yaocirc = std::make_shared<CircuitWrapper>(sharings[S_YAO]->GetCircuitBuildRoutine());
    CircuitW_p boolcirc = std::make_shared<CircuitWrapper>(sharings[S_BOOL]->GetCircuitBuildRoutine());

    auto sorted_cycles = readCyclesFromFile(boolcirc, role, bitlen, n_cycles, cycle_length);

    auto largest_set = findSolution(sorted_cycles, boolcirc, yaocirc, bitlen, n_pairs, n_cycles, cycle_length);

    auto result_weight = std::get<0>(largest_set);
    auto result_cycles = std::get<1>(largest_set);

    result_weight = yaocirc->PutOUTGate(result_weight, ALL);
    for(size_t i = 0; i < n_cycles; ++i) {
        for(size_t j = 0; j < cycle_length; ++j) {
            result_cycles[i][j] = yaocirc->PutOUTGate(result_cycles[i][j], ALL);
        }
    }

    party->ExecCircuit();

    auto result_weight_plain = result_weight->get_clear_value<uint32_t>();

    std::vector<std::vector<uint32_t>> result_set_plain;
    for(size_t i = 0; i < n_cycles; ++i) {
        std::vector<uint32_t> cycle;
        for(size_t j = 0; j < cycle_length; ++j) {
            cycle.push_back(result_cycles[i][j]->get_clear_value<uint32_t>());
        }
        result_set_plain.push_back(cycle);
    }

    std::cout << "\n\n#####\t#####\t#####\t#####\t#####" << std::endl;
    std::cout << "Result:\t" << result_weight_plain << std::endl;
    for(auto &cycle : result_set_plain) {
        std::cout << "Cycle\t";
        for(auto vertex : cycle) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }
}


cmp_graph readCompGraphFromFile1(CircuitW_p &ac, e_role role, uint32_t bitlen, uint32_t n_pairs)
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
 

share_p compare_cycles(std::vector<share_p> &cycleA, std::vector<share_p> &cycleB, CircuitW_p &yc, 
                uint32_t bitlen, uint32_t cycle_length) 
{
    share_p s_length = yc->PutCONSGate(cycle_length, bitlen);
    share_p s_zero = yc->PutCONSGate((uint32_t) 0, bitlen);
    share_p s_one = yc->PutCONSGate((uint32_t) 1, bitlen);

    share_p count = yc->PutINGate((uint32_t) 0, bitlen, SERVER);
    
    for(size_t i = 0; i < cycle_length; ++i) {
        for(size_t j = 0; j < cycle_length; ++j) {
            share_p sel = yc->PutEQGate(cycleA[i], cycleB[j]);
            share_p add_one = yc->PutADDGate(count, s_one);
            share_p add_zero = yc->PutADDGate(count, s_zero);
            count = yc->PutMUXGate(add_one, add_zero, sel);
        }
    }

    // At most cycle_length vertices can be the same, as each vertex appears at most once within the same cycle
    return yc->PutEQGate(s_length, count);
}


share_p containsCycle(cycles &all_cycles, std::vector<share_p> &cycle, share_p &weight, CircuitW_p &yc, 
                uint32_t bitlen, uint32_t cycle_length)
{
    share_p isOld = yc->PutINGate((uint32_t) 0, 1, SERVER); 

    for(size_t i = 0; i < all_cycles.size(); ++i) {
        auto current_cycle = std::get<1>(all_cycles[i]);
        
        // Search through all vertices for possible duplicates
        share_p same_vertices = compare_cycles(current_cycle, cycle, yc, bitlen, cycle_length);
        
        // Weights should be different, if vertices are the same
        auto same_weight = yc->PutEQGate(std::get<0>(all_cycles[i]), weight);

        auto new_cycle = yc->PutANDGate(same_vertices, same_weight);
        isOld = yc->PutORGate(isOld, new_cycle);
    }
    return isOld;
}


void findCycles(cycles &all_cycles, cmp_graph &s_cmp_graph, std::vector<uint32_t> &current_cycle, 
                uint32_t current_length, share_p &weight, share_p &valid, CircuitW_p &yc, 
                uint32_t bitlen, uint32_t n_pairs, uint32_t cycle_length)
{   
    share_p s_one = yc->PutCONSGate((uint32_t) 1, bitlen);
    share_p s_zero = yc->PutCONSGate((uint32_t) 0, bitlen);
    share_p s_max = yc->PutCONSGate((uint32_t) n_pairs, bitlen);

    // Base case
    if(current_length == cycle_length) {
        // Add edge from last to first vertex
        weight = yc->PutADDGate(weight, s_cmp_graph[current_cycle.back()][current_cycle[0]]);
        auto sel = yc->PutGTGate(s_cmp_graph[current_cycle.back()][current_cycle[0]], s_zero);
        auto add_one = yc->PutADDGate(valid, s_one);
        auto add_zero = yc->PutADDGate(valid, s_zero);
        valid = yc->PutMUXGate(add_one, add_zero, sel);

        // Check whether or not the current cycle is valid (all vertices connected)
        share_p s_cycle_length = yc->PutINGate(cycle_length, bitlen, SERVER);
        auto valid_cycle = yc->PutEQGate(valid, s_cycle_length);
        
        // Secret share current cycle
        std::vector<share_p> cycle;
         for(size_t l = 0; l < cycle_length; ++l) {
            cycle.push_back(yc->PutINGate(current_cycle[l], bitlen, SERVER));
        }

        // Set weight true weight iff valid is true and isNew is true 
        auto cycle_weight = yc->PutMUXGate(weight, s_zero, valid_cycle);

        // Add cycle with correspondig weight
        all_cycles.push_back(std::make_tuple(cycle_weight, cycle));

        // Remove last edge
        auto sub_one = yc->PutSUBGate(valid, s_one);
        auto sub_zero = yc->PutSUBGate(valid, s_zero);
        valid = yc->PutMUXGate(sub_one, sub_zero, sel);
        weight = yc->PutSUBGate(weight, s_cmp_graph[current_cycle.back()][current_cycle[0]]);

    } else {    // Go deeper
        for(size_t i = 0; i < n_pairs; ++i) {
            if(current_length == 0) { // Start of a cycle
                // Add first vertex
                current_cycle.push_back(i);

                // Set weight and valid to zero
                weight = yc->PutINGate((uint32_t) 0, bitlen, SERVER);
                valid = yc->PutINGate((uint32_t) 0, bitlen, SERVER);

                // Recursive step with updated values
                findCycles(all_cycles, s_cmp_graph, current_cycle, current_length + 1, weight, valid, 
                            yc, bitlen, n_pairs, cycle_length);

                // Remove vertex
                current_cycle.pop_back();
            } else {
                // For even cycle lengths greater than 3, we exclude duplicates (e.g., 1->2->1->2) which were count in 
                // compute number of cycles. By doing so, the actual number of cycles with the given length is smaller than 
                // the number of cycles found. This results in better memory usage during the recursion, however worse memory usage 
                // during the construction of the maximum set
                if(std::find(current_cycle.begin(), current_cycle.end(), i) != current_cycle.end()) {
                    continue;
                }
                // Add new edge 
                weight = yc->PutADDGate(weight, s_cmp_graph[current_cycle.back()][i]);
                // Check whether the current last vertex and the new one are connected
                auto sel = yc->PutGTGate(s_cmp_graph[current_cycle.back()][i], s_zero);
                auto add_one = yc->PutADDGate(valid, s_one);
                auto add_zero = yc->PutADDGate(valid, s_zero);
                valid = yc->PutMUXGate(add_one, add_zero, sel);
                // Add new vertex
                current_cycle.push_back(i);

                // Recursive step with updated values
                findCycles(all_cycles, s_cmp_graph, current_cycle, current_length + 1, weight, valid, 
                            yc, bitlen, n_pairs, cycle_length);

                // Remove edge and vertex
                current_cycle.pop_back();
                auto sub_one = yc->PutSUBGate(valid, s_one);
                auto sub_zero = yc->PutSUBGate(valid, s_zero);
                valid = yc->PutMUXGate(sub_one, sub_zero, sel);
                weight = yc->PutSUBGate(weight, s_cmp_graph[current_cycle.back()][i]);
            }
        }
    }
}


cycles kNNSort(cycles &all_cycles, CircuitW_p &yc, uint32_t bitlen, uint32_t n_pairs, uint32_t n_all_cycles, 
                uint32_t n_cycles, uint32_t cycle_length)
{                
    std::vector<share_p> sorted_weights;
    cycles_set sorted_cycles;

    for(size_t i = 0; i < n_cycles + 1; ++i) {
        sorted_weights.push_back(yc->PutINGate((uint32_t) 0, bitlen, SERVER));
        std::vector<share_p> vertices;
        for(size_t j = 0; j < cycle_length; ++j) {
            vertices.push_back(yc->PutCONSGate(n_pairs, bitlen));
        }
        sorted_cycles.push_back(vertices);
    }

    for(size_t i = 0; i < n_all_cycles; ++i) {
        sorted_weights[n_cycles] = std::get<0>(all_cycles[i]);
        sorted_cycles[n_cycles] = std::get<1>(all_cycles[i]);

        for(size_t j = n_cycles; j > 0; --j) {
            auto gt = yc->PutGTGate(sorted_weights[j], sorted_weights[j-1]);
            
            // Cond swap, higher value is at position 1 of res array
            auto res = yc->PutCondSwapGate(sorted_weights[j], sorted_weights[j-1], gt);
            sorted_weights[j-1] = res[1];
            sorted_weights[j] = res[0];
            for(size_t k = 0; k < cycle_length; ++k) {
                res = yc->PutCondSwapGate(sorted_cycles[j][k], sorted_cycles[j-1][k], gt);
                sorted_cycles[j-1][k] = res[1];
                sorted_cycles[j][k] = res[0];
            }
        }
    }

    // Combine sorted cycles and weights
    cycles result;
    for(size_t i = 0; i < n_cycles; ++i) {
        result.push_back(std::make_tuple(sorted_weights[i], sorted_cycles[i]));
    } 
    return result;
}


cycles removeDuplicates(cycles &sorted_cycles, CircuitW_p &yc, uint32_t bitlen, uint32_t n_pairs, 
                    uint32_t n_cycles, uint32_t cycle_length)
{  
    auto s_one = yc->PutCONSGate((uint32_t) 1, bitlen);
    auto s_zero = yc->PutCONSGate((uint32_t) 0, bitlen);
    cycles cycles_duplicates;

    for(size_t i = 0; i < n_cycles; ++i) {
        auto cycle = std::get<1>(sorted_cycles[i]);
        auto weight = std::get<0>(sorted_cycles[i]);
        auto isDuplicate = yc->PutINGate((uint32_t) 0, bitlen, SERVER);
        for(size_t j = 0; j < i; ++j) {
            auto current_cycle = std::get<1>(cycles_duplicates[j]);

            // Idea: A cycle exists only once with a specific order, so we shift the cycle to find the same cycle but with different starting vertices
            // To find all possible shifts, we have to check for all shifts from 1 to cycle_length - 1
            for(size_t k = 1; k < cycle_length; ++k) {
                auto duplicate = yc->PutINGate((uint32_t) 1, bitlen, SERVER);
                for(size_t l = 0; l < cycle_length; ++l) {
                    auto same = yc->PutEQGate(cycle[l], current_cycle[(k+l) % cycle_length]);
                    duplicate = yc->PutANDGate(duplicate, same);
                }
                isDuplicate = yc->PutORGate(isDuplicate, duplicate);
            }
        }

        // Use only the valid correct weight for non-duplicates
        auto cycle_weight = yc->PutMUXGate(s_zero, weight, isDuplicate);
        cycles_duplicates.push_back(std::make_tuple(cycle_weight, cycle));
    }

    return kNNSort(cycles_duplicates, yc, bitlen, n_pairs, n_cycles, (int) n_cycles / (int) cycle_length, cycle_length);
}


void prepareWritingCyclesToFile(cycles &sorted_cycles, CircuitW_p &bc, uint32_t n_cycles, uint32_t cycle_length)
{
    for(size_t i = 0; i < n_cycles; ++i) {
        auto weight = std::get<0>(sorted_cycles[i]);
        weight = bc->PutSharedOUTGate(weight);
        auto cycle = std::get<1>(sorted_cycles[i]);
        for(size_t j = 0; j < cycle_length; ++j) {
            cycle[j] = bc->PutSharedOUTGate(cycle[j]);
        }
    }    
}


void writeCyclesToFile(cycles &sorted_cycles, e_role role, uint32_t n_cycles, uint32_t cycle_length)
{
    std::string roleString = "Client";
    if(role == SERVER) {
        roleString = "Server";
    }
    std::ofstream outfile("../data/output/cycles_"+roleString+".json");

    json jycles;
    jycles["n_cycles"] = n_cycles;
    jycles["cycle_length"] = cycle_length;
    json jycle_array = json::array();

    for(size_t i = 0; i < n_cycles; ++i) {
        json jycle;
        auto key = std::to_string(i);
        auto weight = std::get<0>(sorted_cycles[i])->get_clear_value<uint32_t>();
        jycle["weight"] =  weight;
        auto cycle = std::get<1>(sorted_cycles[i]);
        for(size_t j = 0; j < cycle_length; ++j) {
            jycle["vertices"].push_back(cycle[j]->get_clear_value<uint32_t>());
        }
        jycle_array.push_back(jycle);
    }

    jycles["cycles"] = jycle_array;

    outfile << jycles;
}


cycles readCyclesFromFile(CircuitW_p &bc, e_role role, uint32_t bitlen, uint32_t n_cycles, uint32_t cycle_length)
{
    cycles s_cycles;

    std::string roleString = "Client";
    if(role == SERVER) {
        roleString = "Server";
    }
    std::ifstream infile("../data/output/cycles_"+roleString+".json");

    json jycles;
    infile >> jycles;

    if( ((uint32_t) jycles["n_cycles"]) != n_cycles || ((uint32_t) jycles["cycle_length"]) != cycle_length) {
        std::cout << "Number of cycles or cycle length do not match!" <<std::endl;
        return s_cycles;
    }

    // Read cycles
    for(auto &jycle: jycles["cycles"]) {
        auto weight = bc->PutSharedINGate((uint32_t) jycle["weight"], bitlen);
        auto vertices = jycle["vertices"];
        std::vector<share_p> cycle;
        for(size_t i = 0; i < cycle_length; ++i) {
            cycle.push_back(bc->PutSharedINGate((uint32_t) vertices[i], bitlen));
        }
        s_cycles.push_back(std::make_tuple(weight, cycle));
    }

    return s_cycles;
}


share_p disjoint_set(cycles_set set_cycles, std::vector<share_p> cycle, CircuitW_p bc, uint32_t bitlen, 
                uint32_t cycle_count, uint32_t cycle_length)
{
    std::vector<share_p> notD;
    for(size_t i = 0; i < cycle_count; ++i) {
        auto current_cycle = set_cycles[i];
        for(size_t j = 0; j < cycle_length; ++j) {
            for(size_t k = 0; k < cycle_length; ++k) {
                notD.push_back(bc->PutEQGate(current_cycle[j], cycle[k]));
            }
        }
        // Using the intermediary results in notD, we compute the whether two cycles are disjoint
        while(notD.size() > 1) {
            uint32_t i = 0;
            for(size_t j = 0; j < notD.size();) {
                if(j + 1 >= notD.size()) {
                    notD[i++] = notD[j++];
                } else {
                    notD[i++] = bc->PutORGate(notD[j], notD[j+1]);
                }
                j += 2;
            }
            notD.resize(i);
        }
    }
    return bc->PutINVGate(notD[0]);
}


result_set findMaximumSet(std::vector<cycles_set> &all_cycle_sets, std::vector<share_p> &set_weights, 
                CircuitW_p &yc, uint32_t bitlen, uint32_t n_pairs, uint32_t n_cycles, 
                uint32_t cycle_length)
{
    std::vector<cycles_set> sets;
    std::vector<share_p> weights;

    // We only want the set with the highest weight
    for(size_t i = 0; i < 2; ++i) {
        weights.push_back(yc->PutINGate((uint32_t) 0, bitlen, SERVER));
        cycles_set cycles;
        for(size_t j = 0; j < n_cycles; ++j) {
            std::vector<share_p> vertices;
            for(size_t k = 0; k < cycle_length; ++k) {
                vertices.push_back(yc->PutINGate(n_pairs, bitlen, SERVER));
            }
            cycles.push_back(vertices);
        }
        sets.push_back(cycles);
    }

    for(size_t i = 0; i < n_cycles; ++i) {
        weights[1] = set_weights[i];
        sets[1] = all_cycle_sets[i];
        
        auto gt = yc->PutGTGate(weights[1], weights[0]);
            
        // Cond swap, higher value is at position 1 of res array
        auto res = yc->PutCondSwapGate(weights[1], weights[0], gt);
        weights[0] = res[1];
        weights[1] = res[0];
        // Swap each vertex of each cycle
        for(size_t j = 0; j < n_cycles; ++j) {
            for(size_t k = 0; k < cycle_length; ++k) {
                res = yc->PutCondSwapGate(sets[1][j][k], sets[0][j][k], gt);
                sets[0][j][k] = res[1];
                sets[1][j][k] = res[0];
            }
        }
    }

    return result_set(weights[0], sets[0]); 
}


result_set findSolution(cycles &all_cycles, CircuitW_p &bc, CircuitW_p &yc, uint32_t bitlen, uint32_t n_pairs,
                uint32_t n_cycles, uint32_t cycle_length)
{   
    share_p zero = bc->PutCONSGate((uint32_t) 0, bitlen);
    share_p max_pair = bc->PutCONSGate(n_pairs, bitlen);
    std::vector<cycles_set> all_cycle_sets;
    std::vector<share_p> set_weights;

    for(size_t i = 0; i < n_cycles; ++i) {
        cycles_set current_set;
        current_set.push_back(std::get<1>(all_cycles[i]));
        auto weight = std::get<0>(all_cycles[i]);

        // Keep track of the number cycles in current_set
        auto cycle_count = 1;
        for(size_t j = 0; j < n_cycles; ++j) {
            if(i == j) {
                continue;
            }
          
            auto current_cycle = std::get<1>(all_cycles[j]);
            share_p disjoint = disjoint_set(current_set, current_cycle, bc, bitlen, cycle_count, cycle_length);

            std::vector<share_p> vertices;
            for(size_t k = 0; k < cycle_length; ++k) {
                vertices.push_back(bc->PutMUXGate(current_cycle[k], max_pair, disjoint));
            }
            current_set.push_back(vertices);
            auto add = bc->PutADDGate(weight, std::get<0>(all_cycles[j]));
            auto dont_add = bc->PutADDGate(weight, zero);
            weight = bc->PutMUXGate(add, dont_add, disjoint);
            cycle_count++;
        }
    
        all_cycle_sets.push_back(current_set);
        set_weights.push_back(weight); 
    }

    // Convert to yao sharing
    for(size_t i = 0; i < n_cycles; ++i) {
        set_weights[i] = yc->PutB2YGate(set_weights[i]);
        for(size_t j = 0; j < n_cycles; ++j) {
            for(size_t k = 0; k < cycle_length; ++k) {
                all_cycle_sets[i][j][k] = yc->PutB2YGate(all_cycle_sets[i][j][k]);
            }
        }
    }

    // Return the set with the maximum weight
    return findMaximumSet(all_cycle_sets, set_weights, yc, bitlen, n_pairs, n_cycles, cycle_length);
}



// Plain Implementation
void plain_solution(uint32_t n_pairs, uint32_t cycle_length, uint32_t n_all_cycles, uint32_t n_cycles)
{   
    // std::cout << "Start plain solution" << std::endl;

    auto graph_plain = readCompGraphFromFilePlain1(n_pairs);

    // std::cout << "Read graph." << std::endl;

    std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> all_cycles;
    std::vector<uint32_t> current_cycle;
    uint32_t weight = 0;
    uint32_t valid = 0;
    uint32_t current_length = 0;
    all_cycles = findCyclesPlain(all_cycles, graph_plain, current_cycle, current_length, weight, 
            valid, n_pairs, cycle_length); 

    // std::cout << "Found cycles-> size\t" << all_cycles.size() << std::endl;

    all_cycles = sortCyclesPlain(all_cycles, n_all_cycles);

    auto result = findSolutionPlain(all_cycles, n_pairs, n_cycles, cycle_length);

    std::cout << "Expected:\t" << std::get<0>(result) << std::endl;
    for(auto &cycle : std::get<1>(result)) {
        std::cout << "Cycle\t";
        for(auto vertex : cycle) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }
}


std::vector<std::vector<uint32_t>> readCompGraphFromFilePlain1(uint32_t n_pairs) 
{    
    std::vector<std::vector<uint32_t>> comp_graph;

    std::ifstream infile("../data/output/comp_graph_plain.json");
    std::string tmp;

    json jomp_graph;
    infile >> jomp_graph;

    if(jomp_graph["n_pairs"] != n_pairs) {
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


std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> sortCyclesPlain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> cycles,
            uint32_t n_all_cycles) 
{
    int i, j;
    for(i = 1; i < n_all_cycles; ++i) {
        auto key = cycles[i];

        j = i - 1;

        while(j >= 0 && ((int) std::get<0>(cycles[j]) < (int) std::get<0>(key))) {
            cycles[j+1] = cycles[j];
            j = j - 1;
        }
        cycles[j+1] = key;
    }
    return cycles;
}


bool containsCyclePlain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> &all_cycles, std::vector<uint32_t> &cycle, 
            uint32_t weight, uint32_t cycle_length)
{
    bool isOld = false;
    
    for(size_t i = 0; i < all_cycles.size(); ++i) {
        uint32_t count = 0;
        auto current_cycle = std::get<1>(all_cycles[i]);
        for(size_t j = 1; j < cycle_length; ++j) {
            for(size_t k = 0; k < cycle_length; ++k) {
                count = current_cycle[(k+j) % cycle_length] == cycle[k] ? count+1 : count;
            }
        }
        isOld |= count == cycle_length;
    }

    return isOld;
}


std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> findCyclesPlain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> all_cycles,
            std::vector<std::vector<uint32_t>> cmp_graph, std::vector<uint32_t> current_cycle, uint32_t current_length, uint32_t weight, 
            uint32_t valid, uint32_t n_pairs, uint32_t cycle_length) 
{
    // Base case
    if(current_length == cycle_length) {
        // Add edge from last to first vertex
        weight += cmp_graph[current_cycle.back()][current_cycle[0]];
        valid = cmp_graph[current_cycle.back()][current_cycle[0]] > 0 ? valid + 1 : valid;

        // Check if cycles is new
        auto isNew = !containsCyclePlain(all_cycles, current_cycle, weight, cycle_length);

        std::vector<uint32_t> out_cycle;
        for(size_t i = 0; i < cycle_length; ++i) {
            out_cycle.push_back((valid == cycle_length && isNew)? current_cycle[i] : n_pairs);
        }
        // Set cycle weight
        auto cycle_weight = (valid == cycle_length && isNew)? weight : 0;
        
        all_cycles.push_back(std::make_tuple(cycle_weight, out_cycle));

        // Remove last edge
        weight -= cmp_graph[current_cycle.back()][current_cycle[0]];
        valid = cmp_graph[current_cycle.back()][current_cycle[0]] > 0 ? valid - 1 : valid;
    } else {    // Go deeper
        for(size_t i = 0; i < n_pairs; ++i) {
            if(current_length == 0) {
                // Add first vertex
                current_cycle.push_back(i);

                // Set weight and valid to zero
                weight = 0;
                valid = 0;

                // Recursive step with updated values
                all_cycles = findCyclesPlain(all_cycles, cmp_graph, current_cycle, current_length + 1,
                                        weight, valid, n_pairs, cycle_length);
            
                // Remove vertex
                current_cycle.pop_back();
            } else {
                // If vertex is already in the current cycle, continue
                if(std::find(current_cycle.begin(), current_cycle.end(), i) != current_cycle.end()) {
                    continue;
                }
                // Add new edge 
                weight += cmp_graph[current_cycle.back()][i];
                valid = cmp_graph[current_cycle.back()][i] > 0 ? valid + 1 : valid;
                
                // Add new vertex
                current_cycle.push_back(i);

                // Recursive step with updated values
                all_cycles = findCyclesPlain(all_cycles, cmp_graph, current_cycle, current_length + 1,
                                        weight, valid, n_pairs, cycle_length);

                // Remove edge and vertex
                current_cycle.pop_back();
                weight -= cmp_graph[current_cycle.back()][i];
                valid = cmp_graph[current_cycle.back()][i] > 0 ? valid - 1 : valid;
            }
        }
    }
    return all_cycles;
}


bool disjoint_set_plain(std::vector<std::vector<uint32_t>> cycle_set, std::vector<uint32_t> cycle, 
            uint32_t cycle_count, uint32_t cycle_length) {
    
    bool not_disjoint = false;
    for(size_t i = 0; i < cycle_count; ++i) {
        auto current_cycle = cycle_set[i];
        for(size_t j = 0; j < cycle_length; ++j) {
            for(size_t k = 0; k < cycle_length; ++k) {
                not_disjoint = not_disjoint || current_cycle[j] == cycle[k];   
            }
        }
    }
    return !not_disjoint;
}


std::tuple<uint32_t, std::vector<std::vector<uint32_t>>> find_maximum_set_plain(std::vector<std::vector<std::vector<uint32_t>>> cycle_sets, 
        std::vector<uint32_t> set_weight, uint32_t n_cycles) {
    
    uint32_t max = 0;
    for(size_t i = 1; i < n_cycles; ++i) {
        if (set_weight[i] > set_weight[max]) {
            max = i;
        }
    }

    auto tmp_set = cycle_sets[0];
    cycle_sets[0] = cycle_sets[max];
    cycle_sets[max] = tmp_set;
    auto tmp = set_weight[0];
    set_weight[0] = set_weight[max];
    set_weight[max] = tmp;

    return std::make_tuple(set_weight[0], cycle_sets[0]);
}


std::tuple<uint32_t, std::vector<std::vector<uint32_t>>> findSolutionPlain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> cycles, 
        uint32_t n_pairs, uint32_t n_cycles, uint32_t cycle_length) {

    std::vector<std::vector<std::vector<uint32_t>>> cycle_sets;
    std::vector<uint32_t> set_weight;
    std::vector<uint32_t> dummy_cycle;
    for(size_t i = 0; i < cycle_length; ++i) {
        dummy_cycle.push_back(n_pairs);
    }

    for(size_t i = 0; i < n_cycles; ++i) {
        std::vector<std::vector<uint32_t>> current_set;
        current_set.push_back(std::get<1>(cycles[i]));
        uint32_t weight = std::get<0>(cycles[i]);
        auto cycle_count = 1;
        for(size_t j = 0; j < n_cycles; ++j) {
            if(i == j) {
                continue;
            }
            bool disjoint = disjoint_set_plain(current_set, std::get<1>(cycles[j]), cycle_count, cycle_length);

            current_set.push_back(disjoint ? std::get<1>(cycles[j]) : dummy_cycle);
            weight += disjoint ? std::get<0>(cycles[j]) : 0;
            cycle_count++;
        }
        cycle_sets.push_back(current_set);
        set_weight.push_back(weight);
    }

    // Print sets
    // for(size_t i = 0; i < n_cycles; ++i) {
    //     // auto cycles = remove_duplicates(cycle_sets[i], n_cycles, cycle_length);
    //     if(cycles.size() > 1) {
    //         std::cout << "##### New set #####\n" << std::endl;
    //         std::cout << "Weight\t" << set_weight[i] << std::endl;
    //         for(auto &cycle : cycles) {
    //             std::cout << "# New cycle #" << std::endl;
    //             for(auto &vertex : cycle) {
    //                 std::cout << "Vertex:\t" << vertex << std::endl;
    //             }
    //         }
    //         std::cout << std::endl;
    //     }
    // // }

    auto result = find_maximum_set_plain(cycle_sets, set_weight, n_cycles);
    return result;
}
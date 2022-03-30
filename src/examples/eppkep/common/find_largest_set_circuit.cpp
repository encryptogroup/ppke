#include "find_largest_set_circuit.hpp"
#include "../../../abycore/circuit/booleancircuits.h"
#include "../../../abycore/sharing/sharing.h"

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

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
    
    solve(role, address, port, seclvl, bitlen, nthreads, mt_alg, n_pairs, cycle_length, n_all_cycles,
                n_cycles);
    
    std::cout << "\n#####\t#####\t#####\t#####\t#####" << std::endl;
    std::cout << "#####\t#####\t#####\t#####\t#####\n" << std::endl;

    plain_solution(n_pairs, cycle_length, n_all_cycles, n_cycles); 
}   


void solve(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t bitlen, 
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

    cmp_graph s_cmp_graph_gmw;
    // Convert sharing to boolean sharing A2B
    for(size_t i = 0; i < n_pairs; ++i) {
        std::vector<share_p> current_row;
        for(size_t j = 0; j < n_pairs; ++j) {
            current_row.push_back(boolcirc->PutA2BGate(s_cmp_graph[i][j], yaocirc->circ_));
        }
        s_cmp_graph_gmw.push_back(current_row);
    }

    // Prepare for find cycles
    cycles all_cycles;
    share_p weight = arithmcirc->PutINGate((uint32_t) 0, bitlen, SERVER);
    share_p valid = boolcirc->PutINGate((uint32_t) 0, bitlen, SERVER);
    std::vector<uint32_t> current_cycle;

    all_cycles = find_cycles(all_cycles, s_cmp_graph, s_cmp_graph_gmw, current_cycle, 0, weight, 
                                valid, arithmcirc, boolcirc, yaocirc, bitlen, n_pairs, cycle_length);

    cycles all_cycles_yao;
    for(size_t i = 0; i < n_all_cycles; ++i) {
        std::vector<share_p> current_row;
        auto cycle = std::get<1>(all_cycles[i]);
        for(size_t j = 0; j < cycle_length; ++j) {
            current_row.push_back(yaocirc->PutB2YGate(cycle[j]));
        }
        auto weight = yaocirc->PutB2YGate(std::get<0>(all_cycles[i]));
        all_cycles_yao.push_back(std::make_tuple(weight, current_row));
    }

    auto sorted_cycles = kNNSort(all_cycles_yao, yaocirc, bitlen, n_pairs, n_all_cycles, n_cycles, cycle_length);

    cycles sorted_cycles_gmw;
    for(size_t i = 0; i < n_cycles; ++i) {
        std::vector<share_p> current_row;
        auto cycle = std::get<1>(sorted_cycles[i]);
        for(size_t j = 0; j < cycle_length; ++j) {
            current_row.push_back(boolcirc->PutY2BGate(cycle[j]));
        }
        auto weight = boolcirc->PutY2BGate(std::get<0>(sorted_cycles[i]));
        sorted_cycles_gmw.push_back(std::make_tuple(weight, current_row));
    }

    auto result_set = find_solution(sorted_cycles_gmw, boolcirc, yaocirc, bitlen, n_cycles, cycle_length);

    auto result_weight = std::get<0>(result_set);
    auto result_cycles = std::get<1>(result_set);

    result_weight = yaocirc->PutOUTGate(result_weight, ALL);
    for(size_t i = 0; i < n_cycles; ++i) {
        for(size_t j = 0; j < cycle_length; ++j) {
            result_cycles[i][j] = yaocirc->PutOUTGate(result_cycles[i][j], ALL);
        }
    }

    party->ExecCircuit();

    // std::cout << "After circuit execution" << std::endl;

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


cmp_graph readCompGraphFromFile1(CircuitW_p ac, e_role role, uint32_t bitlen, uint32_t n_pairs)
{
    cmp_graph s_cmp_graph;

    std::string roleString = "Client";
    if(role == SERVER) {
        roleString = "Server";
    }
    std::ifstream infile("comp_graph_"+roleString+".txt");
    std::string tmp;

    // First,  informational line
    infile >> tmp;
    for(size_t i = 0; i < n_pairs; ++i) {
        size_t count_entries = 0;
        // n_pairs lines, each line contains n_pairs entries
        std::vector<share_p> current_row;
        std::string line;
        infile >> line;
        std::istringstream s_stream(line);
        while(std::getline(s_stream, tmp, ',') && count_entries++ < n_pairs) {
            if(tmp.empty()) {
                continue;
            }
            current_row.push_back(ac->PutSharedINGate(static_cast<uint32_t>(std::stoul(tmp)), bitlen));
        } 
        s_cmp_graph.push_back(current_row);
    }

    return s_cmp_graph;
}
 

cycles find_cycles(cycles all_cycles, cmp_graph s_cmp_graph, cmp_graph s_cmp_graph_gmw, 
                std::vector<uint32_t> current_cycle, uint32_t current_length, share_p weight, share_p valid,
                CircuitW_p ac, CircuitW_p bc, CircuitW_p yc, uint32_t bitlen, uint32_t n_pairs, 
                uint32_t cycle_length)
{   
    share_p s_one = bc->PutCONSGate((uint32_t) 1, bitlen);
    share_p s_zero = bc->PutCONSGate((uint32_t) 0, bitlen);

    // Base case
    if(current_length == cycle_length) {
        // Add edge from last to first vertex
        weight = ac->PutADDGate(weight, s_cmp_graph[current_cycle.back()][current_cycle[0]]);
        auto sel = bc->PutGTGate(s_cmp_graph_gmw[current_cycle.back()][current_cycle[0]], s_zero);
        auto add_one = bc->PutADDGate(valid, s_one);
        auto add_zero = bc->PutADDGate(valid, s_zero);
        valid = bc->PutMUXGate(add_one, add_zero, sel);

        // Check whether or not the current cycle is valid (all vertices connected)
        share_p s_cycle_length = bc->PutINGate(cycle_length, bitlen, SERVER);
        auto add_cycle = bc->PutEQGate(valid, s_cycle_length);
        auto cycle_weight = bc->PutMUXGate(bc->PutA2BGate(weight, yc->circ_), s_zero, add_cycle);

        // Add current cycle to cycles vector with the respective weight
        std::vector<share_p> s_current_cycle;
        for(size_t l = 0; l < cycle_length; ++l) {
            s_current_cycle.push_back(bc->PutINGate(current_cycle[l], bitlen, SERVER));
        }
        all_cycles.push_back(std::make_tuple(cycle_weight, s_current_cycle));

        // Remove last edge
        auto sub_one = bc->PutSUBGate(valid, s_one);
        auto sub_zero = bc->PutSUBGate(valid, s_zero);
        valid = bc->PutMUXGate(sub_one, sub_zero, sel);
        weight = ac->PutSUBGate(weight, s_cmp_graph[current_cycle.back()][current_cycle[0]]);

    } else {    // Go deeper
        for(size_t i = 0; i < n_pairs; ++i) {
            
            if(current_length == 0) {
                // Add first vertex
                current_cycle.push_back(i);

                // Set weight and valid to zero
                weight = ac->PutINGate((uint32_t) 0, bitlen, SERVER);
                valid = bc->PutINGate((uint32_t) 0, bitlen, SERVER);

                // Recursive step with updated values
                all_cycles = find_cycles(all_cycles, s_cmp_graph, s_cmp_graph_gmw, current_cycle, current_length + 1,
                                        weight, valid, ac, bc, yc, bitlen, n_pairs, cycle_length);
                
                // Remove vertex
                current_cycle.pop_back();
            } else {
                // If vertex is already in the current cycle, and cycle_length < 4 continue
                // If cycle_length is 4 or higher, then the result of part 2 includes cycles as (1, 2, 1, 2)
                // So we cannot exclude duplicates in this case or we will have problems later on
                if(cycle_length < 4 && std::find(current_cycle.begin(), current_cycle.end(), i) != current_cycle.end()) {
                    continue;
                }
                // Add new edge 
                weight = ac->PutADDGate(weight, s_cmp_graph[current_cycle.back()][i]);
                // Check whether the current last vertex and the new one are connected
                auto sel = bc->PutGTGate(s_cmp_graph_gmw[current_cycle.back()][i], s_zero);
                auto add_one = bc->PutADDGate(valid, s_one);
                auto add_zero = bc->PutADDGate(valid, s_zero);
                valid = bc->PutMUXGate(add_one, add_zero, sel);
                // Add new vertex
                current_cycle.push_back(i);

                // Recursive step with updated values
                all_cycles = find_cycles(all_cycles, s_cmp_graph, s_cmp_graph_gmw, current_cycle, current_length + 1,
                                        weight, valid, ac, bc, yc, bitlen, n_pairs, cycle_length);

                // Remove edge and vertex
                current_cycle.pop_back();
                auto sub_one = bc->PutSUBGate(valid, s_one);
                auto sub_zero = bc->PutSUBGate(valid, s_zero);
                valid = bc->PutMUXGate(sub_one, sub_zero, sel);
                weight = ac->PutSUBGate(weight, s_cmp_graph[current_cycle.back()][i]);
            }
        }
    }
    return all_cycles;
}


cycles kNNSort(cycles all_cycles, CircuitW_p yc, uint32_t bitlen, uint32_t n_pairs, uint32_t n_all_cycles, 
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


share_p disjoint_set(cycles_set set_cycles, std::vector<share_p> cycle, CircuitW_p bc, uint32_t bitlen, 
                uint32_t cycle_count, uint32_t cycle_length)
{
    share_p not_disjoint = bc->PutINGate((uint32_t) 0, bitlen, SERVER);
    for(size_t i = 0; i < cycle_count; ++i) {
        auto current_cycle = set_cycles[i];
        for(size_t j = 0; j < cycle_length; ++j) {
            for(size_t k = 0; k < cycle_length; ++k) {
                auto tmp  = bc->PutEQGate(current_cycle[j], cycle[k]);
                not_disjoint = bc->PutORGate(not_disjoint, tmp);
            }
        }
    }
    return bc->PutINVGate(not_disjoint);
}


result find_maximum_set(std::vector<cycles_set> all_cycle_sets, std::vector<share_p> set_weights, 
                CircuitW_p yc, uint32_t bitlen, uint32_t n_cycles, uint32_t cycle_length)
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
                vertices.push_back(yc->PutINGate((uint32_t) 100, bitlen, SERVER));
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

    return result(weights[0], sets[0]); 
}


result find_solution(cycles all_cycles, CircuitW_p bc, CircuitW_p yc, uint32_t bitlen, uint32_t n_cycles, 
                uint32_t cycle_length)
{   
    auto zero = bc->PutCONSGate((uint32_t) 0, bitlen);
    std::vector<cycles_set> all_cycle_sets;
    std::vector<share_p> set_weights;

    for(size_t i = 0; i < n_cycles; ++i) {
        cycles_set current_set;
        current_set.push_back(std::get<1>(all_cycles[i]));
        auto weight = std::get<0>(all_cycles[i]);

        // Used to prevent adding invalid cycles to the current set
        auto dummy_cycle = current_set[0];
        // Keep track of the cycles in current_set
        auto cycle_count = 1;
        for(size_t j = 0; j < n_cycles; ++j) {
            if(i == j) {
                continue;
            }
          
            auto current_cycle = std::get<1>(all_cycles[j]);
            share_p disjoint = disjoint_set(current_set, current_cycle, bc, bitlen, cycle_count, cycle_length);

            std::vector<share_p> vertices;
            for(size_t k = 0; k < cycle_length; ++k) {
                vertices.push_back(bc->PutMUXGate(current_cycle[k], dummy_cycle[k], disjoint));
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
    return find_maximum_set(all_cycle_sets, set_weights, yc, bitlen, n_cycles, cycle_length);
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
    all_cycles = find_cycles_plain(all_cycles, graph_plain, current_cycle, current_length, weight, 
            valid, n_pairs, cycle_length); 

    // std::cout << "Found cycles-> size\t" << all_cycles.size() << std::endl;

    all_cycles = sort_cycles(all_cycles, n_all_cycles);

    // std::cout << "Sorted cycles" << std::endl;

    auto result = find_solution_plain(all_cycles, n_cycles, cycle_length);

    std::cout << "Expected:\t" << std::get<0>(result) << std::endl;
    for(auto &cycle : std::get<1>(result)) {
        std::cout << "Cycle\t";
        for(auto vertex : cycle) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }
}


std::vector<std::vector<uint32_t>> readCompGraphFromFilePlain1(uint32_t n_pairs) {
    
    std::vector<std::vector<uint32_t>> comp_graph;
    std::ifstream infile("comp_graph_plain.txt");
    std::string tmp;

    // First,  informational line
    infile >> tmp;
    for(size_t i = 0; i < n_pairs; ++i) {
        // n_pairs lines, each line contains n_pairs entries
        size_t count_entries = 0;
        std::vector<uint32_t> current_row;
        std::string line;
        infile >> line;
        std::istringstream s_stream(line);
        while(std::getline(s_stream, tmp, ',') && count_entries++ < n_pairs) {
            if(tmp.empty()) {
                continue;
            }
            current_row.push_back(static_cast<uint32_t>(std::stoul(tmp)));
        }
        comp_graph.push_back(current_row);
    }

    return comp_graph;
}


std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> sort_cycles(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> cycles,
            uint32_t n_all_cycles) {

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


std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> find_cycles_plain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> all_cycles,
            std::vector<std::vector<uint32_t>> cmp_graph, std::vector<uint32_t> current_cycle, uint32_t current_length, uint32_t weight, 
            uint32_t valid, uint32_t n_pairs, uint32_t cycle_length) 
{
    // Base case
    if(current_length == cycle_length) {
        // Add edge from last to first vertex
        weight += cmp_graph[current_cycle.back()][current_cycle[0]];
        valid = cmp_graph[current_cycle.back()][current_cycle[0]] > 0 ? valid + 1 : valid;

        // Add current cycle to cycles vector with the respective weight
        auto cycle_weight = valid == cycle_length ? weight : 0;
        all_cycles.push_back(std::make_tuple(cycle_weight, current_cycle));

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
                all_cycles = find_cycles_plain(all_cycles, cmp_graph, current_cycle, current_length + 1,
                                        weight, valid, n_pairs, cycle_length);
            
                // Remove vertex
                current_cycle.pop_back();
            } else {
                // If vertex is already in the current cycle, continue
                if( cycle_length < 4 && std::find(current_cycle.begin(), current_cycle.end(), i) != current_cycle.end()) {
                    continue;
                }
                // Add new edge 
                weight += cmp_graph[current_cycle.back()][i];
                valid = cmp_graph[current_cycle.back()][i] > 0 ? valid + 1 : valid;
                
                // Add new vertex
                current_cycle.push_back(i);

                // Recursive step with updated values
                all_cycles = find_cycles_plain(all_cycles, cmp_graph, current_cycle, current_length + 1,
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


std::tuple<uint32_t, std::vector<std::vector<uint32_t>>> find_solution_plain(std::vector<std::tuple<uint32_t, std::vector<uint32_t>>> cycles, 
        uint32_t n_cycles, uint32_t cycle_length) {

    std::vector<std::vector<std::vector<uint32_t>>> cycle_sets;
    std::vector<uint32_t> set_weight;
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

            current_set.push_back(disjoint ? std::get<1>(cycles[j]) : current_set[0]);
            weight += disjoint ? std::get<0>(cycles[j]) : 0;
            cycle_count++;
        }
        cycle_sets.push_back(current_set);
        set_weight.push_back(weight);
    }

    // Print sets
    // for(size_t i = 0; i < n_cycles; ++i) {
    //     auto cycles = remove_duplicates(cycle_sets[i], n_cycles, cycle_length);
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
    // }

    auto result = find_maximum_set_plain(cycle_sets, set_weight, n_cycles);
    return result;
}
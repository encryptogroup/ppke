#include "matching_circuit.hpp"
#include "../../../abycore/circuit/booleancircuits.h"
#include "../../../abycore/sharing/sharing.h"

#include "patient.hpp"
#include "pair.hpp"
#include "plain_data.hpp"

#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

void matching_circuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, std::string path_to_data, 
        const bool factors, uint32_t n_pairs, uint32_t n_hla) 
{
    BuildMatchingCircuit(role, address, port, seclvl, bitlen, nthreads, mt_alg, path_to_data,
            factors, n_pairs, n_hla);

    testing(path_to_data, factors, n_pairs, n_hla);
}


void BuildMatchingCircuit(e_role role, const std::string& address, uint16_t port, seclvl seclvl,
		        uint32_t bitlen, uint32_t nthreads, e_mt_gen_alg mt_alg, std::string path_to_data, 
                const bool factors, uint32_t n_pairs, uint32_t n_hla)
{
    ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads,
			mt_alg);

	std::vector<Sharing*>& sharings = party->GetSharings();

    // Create circuits
	CircuitW_p arithmcirc = std::make_shared<CircuitWrapper>(sharings[S_ARITH]->GetCircuitBuildRoutine());
    CircuitW_p yaocirc = std::make_shared<CircuitWrapper>(sharings[S_YAO]->GetCircuitBuildRoutine());
    CircuitW_p boolcirc = std::make_shared<CircuitWrapper>(sharings[S_BOOL]->GetCircuitBuildRoutine());

    // Initialisation of constants -> Has to be done secret shared 

    // Initialise constants for HLA match
    share_p s_best_hla = boolcirc->PutCONSGate((uint32_t) 5, bitlen);
    share_p s_good_hla = boolcirc->PutCONSGate((uint32_t) 3, bitlen);
    share_p s_ok_hla = boolcirc->PutCONSGate((uint32_t) 1, bitlen);
    share_p s_zero_hla = boolcirc->PutCONSGate((uint32_t) 0, bitlen);
	
    // Initialise constants for ABO match
    share_p s_best_abo = boolcirc->PutCONSGate((uint32_t) 3, bitlen);
    share_p s_zero_abo = boolcirc->PutCONSGate((uint32_t) 0, bitlen);

    // Initialise constants for age match
    share_p s_best_age = boolcirc->PutCONSGate((uint32_t) 2, bitlen);
    share_p s_good_age = boolcirc->PutCONSGate((uint32_t) 1, bitlen);
    share_p s_zero_age = boolcirc->PutCONSGate((uint32_t) 0, bitlen);

    //Initialise constants for sex match
    share_p s_best_sex = boolcirc->PutCONSGate((uint32_t) 2, bitlen);
    share_p s_good_sex = boolcirc->PutCONSGate((uint32_t) 1, bitlen);
    share_p s_zero_sex = boolcirc->PutCONSGate((uint32_t) 0, bitlen);

    // Initialise constants for size match
    share_p s_good_size = boolcirc->PutCONSGate((uint32_t) 1, bitlen);
    share_p s_zero_size = boolcirc->PutCONSGate((uint32_t) 0, bitlen);

    auto pairs = read_pair_data(boolcirc, role, bitlen, path_to_data, n_pairs);

    // std::cout << "Reading successful" << std::endl;

    // Can be done variable to via command line or similar
    std::vector<uint32_t> weights = {1, 1, 1, 1, 1};
    std::vector<share_p> s_weights;
    for(size_t i = 0; i < weights.size(); ++i) {
        s_weights.push_back(arithmcirc->PutCONSGate(weights[i], bitlen));
    }
    
    // Compute the results of the each biomedical factor
    std::vector<cmp_graph> matching_results;
    auto s_zero = boolcirc->PutCONSGate((uint32_t) 0, bitlen);

    for(size_t i = 0; i < n_pairs; ++i) {
        // i-th pair always uses its donor
        cmp_graph current_row;
        std::shared_ptr<Patient> donor = pairs[i]->getDonor();
        for(size_t j = 0; j < n_pairs; ++j) {
            // j-th pair always uses its recipient
            std::vector<share_p> res;
            if(i == j) {
                res.push_back(s_zero);                                               // crossmatch
                res.push_back(s_zero);                                               // ABO match
                if(factors) {
                    res.push_back(s_zero);                                               // HLA match
                    res.push_back(s_zero);                                               // Age match
                    res.push_back(s_zero);                                               // Sex match
                    res.push_back(s_zero);                                               // Size match
                }
            } else {
                std::shared_ptr<Recipient> recipient = pairs[j]->getRecipient();
                auto hla_d = donor->getHLA();
                auto hla_a_r = recipient->getHLA_A();
                res.push_back(BuildCrossmatchCircuit(hla_d, hla_a_r, boolcirc,      // crossmatch
                                n_hla));
                auto abo_d = donor->getBloodGroup();
                auto abo_r = recipient->getBloodGroup();
                res.push_back(BuildABOmatchCircuit(abo_d, abo_r, s_best_abo,        // ABO match
                                s_zero_abo, boolcirc));
                if(factors) {
                    auto hla_r = recipient->getHLA(); 
                    res.push_back(BuildHLAmatchCircuit(hla_d, hla_r, s_best_hla,        // HLA match   
                                s_good_hla, s_ok_hla, s_zero_hla, boolcirc));
                    auto age_d = donor->getAge();
                    auto age_r = recipient->getAge(); 
                    res.push_back(BuildAgematchCircuit(age_d, age_r, s_best_age,        // Age match
                                s_good_age, s_zero_age, boolcirc));
                    auto sex_d = donor->getSex();
                    auto sex_r = recipient->getSex();
                    res.push_back(BuildSexmatchCircuit(sex_d, sex_r, s_best_sex,        // Sex match
                                s_good_sex, s_zero_sex, boolcirc));
                    auto size_d = donor->getSize();
                    auto size_r = recipient->getSize();             
                    res.push_back(BuildSizematchCircuit(size_d, size_r, s_good_size,    // Size match 
                                s_zero_size, boolcirc)); 
                }
            }
            current_row.push_back(res);
        }
        matching_results.push_back(current_row);
    }

    // std::cout << "Computed results of matching.\n" << std::endl;

    cmp_graph s_comp_graph; // Resulting graph
    for(size_t i = 0; i < n_pairs; ++i) {
        std::vector<share_p> current_row;
        for(size_t j = 0; j < n_pairs; ++j) {
            if(factors) {
                std::vector<share_p> res;
                for(size_t k = 0; k < s_weights.size(); ++k) { // At the 0-th index of matching results 
                                                                // is the result of the cross match stored
                    res.push_back(arithmcirc->PutMULGate(s_weights[k], 
                                            arithmcirc->PutB2AGate(matching_results[i][j][k+1])));
                }

                share_p sum = arithmcirc->PutINGate((uint32_t) 1, bitlen, SERVER);
                share_p zero = arithmcirc->PutINGate((uint32_t) 0, bitlen, SERVER); 
                for(size_t k = 0; k < s_weights.size(); ++k) {
                    sum = arithmcirc->PutADDGate(sum, res[k]);
                }
                // Convert the values into an boolean share for the selection
                share_p sum_bool = boolcirc->PutA2BGate(sum, yaocirc->circ_);
                share_p zero_bool = boolcirc->PutA2BGate(zero, yaocirc->circ_);
                share_p tmp = boolcirc->PutMUXGate(sum_bool, zero_bool, matching_results[i][j][0]);
                // Convert result back to an arithmetic share
                current_row.push_back(arithmcirc->PutB2AGate(tmp));
            } else {
                share_p tmp = boolcirc->PutGTGate(matching_results[i][j][1], s_zero);
                current_row.push_back(arithmcirc->PutB2AGate(boolcirc->PutANDGate(matching_results[i][j][0], tmp)));
            }  
        }
        s_comp_graph.push_back(current_row);
    }
    
    // std::cout << "Computed compatibility graph.\n" << std::endl;

    for(size_t i = 0; i < n_pairs; ++i) {
        for(size_t j = 0; j < n_pairs; ++j) {
            s_comp_graph[i][j] = arithmcirc->PutOUTGate(s_comp_graph[i][j], ALL);
        } 
    }

    prepareWritingCompGraphToFile(s_comp_graph, arithmcirc, n_pairs);

    // std::cout << "Prepared shares for writing to file.\n" << std::endl;

    party->ExecCircuit();

    writeCompGraphToFile(s_comp_graph, role, n_pairs);

    // std::cout << "Wrote compatibility graph to file.\n" << std::endl;

    std::vector<std::vector<uint32_t>> comp_graph;
    for(size_t i = 0; i < n_pairs; ++i) {
        std::vector<uint32_t> current_row;
        for(size_t j = 0; j < n_pairs; ++j) {
            current_row.push_back(s_comp_graph[i][j]->get_clear_value<uint32_t>());
        } 
        comp_graph.push_back(current_row);
    }

    party->Reset();
    delete party;

    std::cout << "\n\n#####\t#####\t#####\t#####\t#####" << std::endl;
    std::cout << "Result:" << std::endl;
    for(size_t i = 0; i < n_pairs; ++i) {
        for(size_t j = 0; j < n_pairs; ++j) {
            std::cout << comp_graph[i][j] << " ";
        }
        std::cout << "\n" << std::endl;
    }
}


std::vector<std::shared_ptr<Pair>> read_pair_data(CircuitW_p bc, e_role role, uint32_t bitlen, 
                std::string &path, uint32_t n_pairs)
{    
    std::string filename = "pair";
	std::vector<std::shared_ptr<Pair>> pairs;

    for(int i = 0; i < n_pairs; ++i) {
        std::ifstream ifs(path + filename + std::to_string(i) + ".txt");
        int n_hla;
        std::vector<std::string> donor_attributes;
        std::vector<std::string> recipient_attributes;
        bool readDonor = true;
        std::string line;

        while(std::getline(ifs, line)) {
            if(line.empty()) {
                continue;
            } else {
                if(line.find("Number of HLA") != std::string::npos) {
                    std::string tmp;
                    std::getline(ifs, tmp);
                    try {
                        n_hla = std::stoi(tmp);
                    } catch(const std::exception &e) {
                        std::cout << "Error when reading the number of hla" << std::endl;
                        return pairs;
                    }
                    continue;
                }

                if(line.find("Recipient") != std::string::npos) {
                    readDonor = false;
                    continue;
                }

                if(readDonor) {
                    if(line.find("Donor") != std::string::npos || line.find("#") != std::string::npos) {
                    continue;
                    } else {
                        donor_attributes.push_back(line);
                    }
                } else {
                    if(line.find("#") != std::string::npos) {
                        continue;
                    } else {
                        recipient_attributes.push_back(line);
                    }
                }
            }

        }

        // Parsing Donor Attributes and creating Donor
        // HLA
        uint32_t* hla_d = new uint32_t[n_hla];
        size_t counter = 0;
        std::istringstream hla_d_stream(donor_attributes[0]); 
        std::string tmp;
        while(std::getline(hla_d_stream, tmp, ' ')) {
            if(tmp.empty()) {
                continue;
            }
            hla_d[counter++] = std::stoi(tmp); 
        }
        
        //ABO
        uint32_t abo_d_tmp = std::stoi(donor_attributes[1]);
        uint32_t* abo_d = new uint32_t[2];
        abo_d[0] = (uint32_t) (abo_d_tmp % 2 == 1);
        abo_d[1] = (uint32_t) (abo_d_tmp >= 2);
        
        //Age
        uint32_t age_d = std::stoi(donor_attributes[2]);
        
        //Sex
        uint32_t sex_d = std::stoi(donor_attributes[3]);
        
        //Size
        uint32_t size_d = std::stoi(donor_attributes[4]);

        //std::cout << "Parsed Donor" << std::endl;

        // Parsing Recipient Attributes
        // HLA
        uint32_t* hla_r = new uint32_t[n_hla];
        counter = 0;
        std::istringstream hla_r_stream(recipient_attributes[0]); 
        while(std::getline(hla_r_stream, tmp, ' ')) {
            if(tmp.empty()) {
                continue;
            }
            hla_r[counter++] = std::stoi(tmp); 
        }

        // HLA_A
        uint32_t* hla_a_r = new uint32_t[n_hla];
        counter = 0;
        std::istringstream hla_a_r_stream(recipient_attributes[1]); 
        while(std::getline(hla_a_r_stream, tmp, ' ')) {
            if(tmp.empty()) {
                continue;
            }
            hla_a_r[counter++] = std::stoi(tmp); 
        }

        //ABO
        uint32_t abo_r_tmp = std::stoi(recipient_attributes[2]);
        uint32_t* abo_r = new uint32_t[2];
        abo_r[0] = (uint32_t) (abo_r_tmp % 2 == 1);
        abo_r[1] = (uint32_t) (abo_r_tmp >= 2);
        
        //Age
        uint32_t age_r = std::stoi(recipient_attributes[3]);
        
        //Sex
        uint32_t sex_r = std::stoi(recipient_attributes[4]);
        
        //Size
        uint32_t size_r = std::stoi(recipient_attributes[5]);

        //std::cout << "Parsed Recipient" << std::endl;
        
        // Create shares for medical data
        share_p s_hla_d, s_hla_r, s_hla_a_r;
        std::vector<share_p> s_abo_d, s_abo_r;
        share_p s_age_d, s_age_r;
        share_p s_sex_d, s_sex_r;
        share_p s_size_d, s_size_r;
        if(role == SERVER) {
            s_hla_d = bc->PutSIMDINGate(n_hla, hla_d, 1, SERVER);
            s_hla_r = bc->PutSIMDINGate(n_hla, hla_r, 1, SERVER);
            s_hla_a_r = bc->PutSIMDINGate(n_hla, hla_a_r, 1, SERVER);

            for(size_t i = 0; i < 2; ++i) {
                s_abo_d.push_back(bc->PutINGate(abo_d[i], 1, SERVER));
                s_abo_r.push_back(bc->PutINGate(abo_r[i], 1, SERVER));
            }

            s_age_d = bc->PutINGate(age_d, 1, SERVER);
            s_age_r = bc->PutINGate(age_r, 1, SERVER);
            s_sex_d = bc->PutINGate(sex_d, 1, SERVER);
            s_sex_r = bc->PutINGate(sex_r, 1, SERVER);
            // bitlen is set to 8 which limits weights to 255KG
            s_size_d = bc->PutINGate(size_d, 8, SERVER);
            s_size_r = bc->PutINGate(size_r, 8, SERVER);
        } else {
            s_hla_d = bc->PutDummySIMDINGate(n_hla, 1);
            s_hla_r = bc->PutDummySIMDINGate(n_hla, 1);
            s_hla_a_r = bc->PutDummySIMDINGate(n_hla, 1);

            for(size_t i = 0; i < 2; ++i) {
                s_abo_d.push_back(bc->PutDummyINGate(1));
                s_abo_r.push_back(bc->PutDummyINGate(1));
            }
           
            s_age_d = bc->PutDummyINGate(1);
            s_age_r = bc->PutDummyINGate(1);
            s_sex_d = bc->PutDummyINGate(1);
            s_sex_r = bc->PutDummyINGate(1);
            s_size_d = bc->PutDummyINGate(8);
            s_size_r = bc->PutDummyINGate(8);
        }

        // Create pairs
        std::shared_ptr<Patient> donor = std::make_shared<Patient>(s_hla_d, s_abo_d, s_age_d, 
                                                                    s_sex_d, s_size_d);
        std::shared_ptr<Recipient> recipient = std::make_shared<Recipient>(s_hla_r, s_hla_a_r, s_abo_r, 
                                                                            s_age_r, s_sex_r, s_size_r);

        pairs.push_back(std::make_shared<Pair>(i, donor, recipient));        
    }

    return pairs;
}


share_p BuildCrossmatchCircuit(share_p s_hla_d, share_p s_hla_a_r, 
        CircuitW_p bc, uint32_t n_hla) 
{    
    // Compute similarities with SIMD shares
	share_p res_and = bc->PutANDGate(s_hla_d, s_hla_a_r);

    // Split simd share for further evaluation
    share_p split = bc->PutSplitterGate(res_and);
    // Initiliase result share with 0
    share_p s_incompatible = bc->PutINGate((uint32_t) 0, 1, SERVER);

    // Reduce the splitted share to a single share
    s_incompatible = bc->PutORGate(share_p(split->get_wire_ids_as_share(0)), s_incompatible);
    for(size_t i = 1; i < n_hla; ++i) {
        s_incompatible = bc->PutORGate(share_p(split->get_wire_ids_as_share(i)), s_incompatible);
    }
    
    // Invert the result of the reduction
    s_incompatible = bc->PutINVGate(s_incompatible);
    
    return s_incompatible;
}


share_p BuildHLAmatchCircuit(share_p s_hla_d, share_p s_hla_r, share_p s_best, share_p s_good, 
        share_p s_ok, share_p s_zero, CircuitW_p bc) 
{            
    // Compute mismatches with simd shares
    share_p mismatches = bc->PutXORGate(s_hla_d, s_hla_r);
    
    // Split share and compute hamming weight distance
    share_p sum = bc->PutSplitterGate(mismatches);
    sum = bc->PutHammingWeightGate(sum);

    // bc->PutPrintValueGate(sum, "HLA sum\t");

    share_p five = bc->PutCONSGate((uint32_t) 5, 32);
    share_p three = bc->PutCONSGate((uint32_t) 3, 32);
    share_p zero = bc->PutCONSGate((uint32_t) 0, 32);
    share_p ok = bc->PutGTGate(five, sum);
    share_p good = bc->PutGTGate(three, sum);
    share_p best = bc->PutEQGate(s_zero, sum);

    //Selection
    share_p sel1 = bc->PutMUXGate(s_ok, s_zero, ok);
    share_p sel2 = bc->PutMUXGate(s_good, sel1, good);
                
    share_p out = bc->PutMUXGate(s_best, sel2, best);
    return out;
}


share_p BuildABOmatchCircuit(std::vector<share_p> &s_abo_d, std::vector<share_p> &s_abo_r, 
        share_p s_best, share_p s_zero, CircuitW_p bc) 
{
    share_p tmp, tmp_l, tmp_r, rLeft, rRight, sel, out;

    // Check for the Same blood group
    tmp_l = bc->PutXORGate(s_abo_d[0], s_abo_r[0]);
    tmp_r = bc->PutXORGate(s_abo_d[1], s_abo_r[1]);
    rLeft = bc->PutORGate(tmp_l, tmp_r);
    rLeft = bc->PutINVGate(rLeft);

    // Reset share_p
    tmp_l.reset();
    tmp_r.reset();
    // Evaluate the right part of the formula
    tmp = bc->PutINVGate(s_abo_d[0]);
    tmp_l = bc->PutANDGate(s_abo_r[1], tmp);
    tmp.reset();
    tmp = bc->PutINVGate(s_abo_d[1]);
    tmp_r = bc->PutANDGate(s_abo_r[0], tmp);
    rRight = bc->PutORGate(tmp_l, tmp_r);
    
    sel = bc->PutORGate(rLeft, rRight);
    out = bc->PutMUXGate(s_best, s_zero, sel);

    // bc->PutPrintValueGate(out, "ABO weight\t");
    return out;
}


share_p BuildAgematchCircuit(share_p s_age_d, share_p s_age_r, share_p s_best, share_p s_good, 
        share_p s_zero, CircuitW_p bc) 
{
    share_p same, tmp, ydor, tmp2, tmp3, out;

    //Determine selection bits
    same = bc->PutEQGate(s_age_d, s_age_r);
    tmp = bc->PutINVGate(s_age_d);
    ydor = bc->PutANDGate(tmp, s_age_r);

    // Compute inner mux gates
    tmp2 = bc->PutMUXGate(s_best, s_good, same);
    tmp3 = bc->PutMUXGate(s_best, s_zero, same);

    //Select correct solution
    out = bc->PutMUXGate(tmp2, tmp3, ydor);

    // bc->PutPrintValueGate(out, "Age weight\t");
    return out;
}


share_p BuildSexmatchCircuit(share_p s_sex_d, share_p s_sex_r, share_p s_best, share_p s_good, 
        share_p s_zero, CircuitW_p bc) 
{
    share_p same, tmp, fdmr, tmp2, tmp3, out;

    //Determine selection bits
    same = bc->PutEQGate(s_sex_d, s_sex_r);
    tmp = bc->PutINVGate(s_sex_r);
    fdmr = bc->PutANDGate(s_sex_d, tmp);

    // Compute inner mux gates
    tmp2 = bc->PutMUXGate(s_best, s_zero, same);
    tmp3 = bc->PutMUXGate(s_best, s_good, same);

    //Select correct solution
    out = bc->PutMUXGate(tmp2, tmp3, fdmr);

    // bc->PutPrintValueGate(out, "Sex weight\t");
    return out;
}


share_p BuildSizematchCircuit(share_p s_size_d, share_p s_size_r, share_p s_good, share_p s_zero, 
        CircuitW_p bc) 
{
    share_p tmp, out;

	tmp = bc->PutGTGate(s_size_r, s_size_d);

    out = bc->PutMUXGate(s_zero, s_good, tmp);
    
    // bc->PutPrintValueGate(out, "Size weight\t");
	return out;
}


void prepareWritingCompGraphToFile(cmp_graph &comp_graph, CircuitW_p ac, uint32_t n_pairs) 
{
    for(size_t i = 0; i < n_pairs; ++i) {
        for(size_t j = 0; j < n_pairs; ++j) {
            comp_graph[i][j] = ac->PutSharedOUTGate(comp_graph[i][j]);
        }
    }
}


void writeCompGraphToFile(cmp_graph comp_graph, e_role role, uint32_t n_pairs) 
{
    std::string roleString = "Client";
    if(role == SERVER) {
        roleString = "Server";
    }
    std::ofstream outfile("comp_graph_"+roleString+".txt");
    outfile << "#CompatibilityGraph" << std::endl;
    for(size_t i = 0; i < n_pairs; ++i) {
        for(size_t j = 0; j < n_pairs; ++j) {
            outfile << comp_graph[i][j]->get_clear_value<uint32_t>() << ",";
        }
        outfile << "\n";
    }
}



// Plain implementation
void testing(std::string path_to_data, const bool factors, uint32_t n_pairs, uint32_t n_hla) 
{   
    auto pairs = read_pair_data_plain(path_to_data, n_pairs);
    std::vector<uint32_t> weights = {1, 1, 1, 1, 1};

    std::vector<std::vector<std::vector<uint32_t>>> matching_results_plain;
    // Computing all matching results
    for(size_t i = 0; i < n_pairs; ++i) {
        // i-th pair always uses its donor
        std::vector<std::vector<uint32_t>> current;
        for(size_t j = 0; j < n_pairs; ++j) {
            // j-th pair always uses its recipient
            std::vector<uint32_t> res;
            if(i == j) {
                res.push_back(0);                                               // crossmatch
                res.push_back(0);                                               // ABO match
                if(factors) {
                    res.push_back(0);                                               // HLA match
                    res.push_back(0);                                               // Age match
                    res.push_back(0);                                               // Sex match
                    res.push_back(0);                                               // Size match
                }
            } else {
                res.push_back(evalCrossmatch(pairs[i]->getHLA_D(), 
                            pairs[j]->getHLA_A_R(), n_hla));                   // crossmatch
                res.push_back(evalABOmatch(pairs[i]->getABO_D(),
                                pairs[j]->getABO_R()));                        // ABO match 
                if(factors) {
                    res.push_back(evalHLAmatch(pairs[i]->getHLA_D(), 
                                pairs[j]->getHLA_R(), n_hla));                     // HLA match
                    res.push_back(evalAgematch(pairs[i]->getAge_D(),
                                pairs[j]->getAge_R()));                            // Age match
                    res.push_back(evalSexmatch(pairs[i]->getSex_D(),
                                pairs[j]->getSex_R()));                              // Sex match
                    res.push_back(evalSizematch(pairs[i]->getSize_D(),
                                pairs[j]->getSize_R()));                           // Size match
                }
            }
            current.push_back(res);
        }
        matching_results_plain.push_back(current);
    }

    std::vector<std::vector<uint32_t>> comp_graph_plain;
    for(size_t i = 0; i < n_pairs; ++i) {
        std::vector<uint32_t> current_row;
        for(size_t j = 0; j < n_pairs; ++j) {
            if(factors) {
                int sum = 1;
                for(size_t k = 0; k < weights.size(); ++k) {
                    sum += weights[k] * matching_results_plain[i][j][k+1];
                }
                current_row.push_back(matching_results_plain[i][j][0] == 1 ? sum : 0);
            } else {
                current_row.push_back(matching_results_plain[i][j][0] > 0 && matching_results_plain[i][j][1] > 0);
            }
        }
        comp_graph_plain.push_back(current_row);
    }

    writeCompGraphToFilePlain(comp_graph_plain, n_pairs);

    std::cout << "#####\t#####\t#####\t#####\t#####" << std::endl;
    std::cout << "Expected:" << std::endl;
    for(size_t i = 0; i < n_pairs; ++i) {
        for(size_t j = 0; j < n_pairs; ++j) {
            std::cout << comp_graph_plain[i][j] << " ";
        }
        std::cout << "\n" << std::endl;
    }
}


std::vector<std::shared_ptr<PlainPair>> read_pair_data_plain(std::string &path, uint32_t n_pairs) 
{
    std::string filename = "pair";
	std::vector<std::shared_ptr<PlainPair>> pairs;

    for(int i = 0; i < n_pairs; ++i) {
        std::ifstream ifs(path + filename + std::to_string(i) + ".txt");
        int n_hla;
        std::vector<std::string> donor_attributes;
        std::vector<std::string> recipient_attributes;
        bool readDonor = true;
        std::string line;

        while(std::getline(ifs, line)) {
            if(line.empty()) {
                continue;
            } else {
                if(line.find("Number of HLA") != std::string::npos) {
                    std::string tmp;
                    std::getline(ifs, tmp);
                    try {
                        n_hla = std::stoi(tmp);
                    } catch(const std::exception &e) {
                        std::cout << "Error when reading the number of hla" << std::endl;
                        return pairs;
                    }
                    continue;
                }

                if(line.find("Recipient") != std::string::npos) {
                    readDonor = false;
                    continue;
                }

                if(readDonor) {
                    if(line.find("Donor") != std::string::npos || line.find("#") != std::string::npos) {
                    continue;
                    } else {
                        donor_attributes.push_back(line);
                    }
                } else {
                    if(line.find("#") != std::string::npos) {
                        continue;
                    } else {
                        recipient_attributes.push_back(line);
                    }
                }
            }

        }

        // Parsing Donor Attributes and creating Donor
        // HLA
        uint32_t* hla_d = new uint32_t[n_hla];
        size_t counter = 0;
        std::istringstream hla_d_stream(donor_attributes[0]); 
        std::string tmp;
        while(std::getline(hla_d_stream, tmp, ' ')) {
            if(tmp.empty()) {
                continue;
            }
            hla_d[counter++] = std::stoi(tmp); 
        }
        
        //ABO
        uint32_t abo_d_tmp = std::stoi(donor_attributes[1]);
        uint32_t* abo_d = new uint32_t[2];
        abo_d[0] = (uint32_t) (abo_d_tmp % 2 == 1);
        abo_d[1] = (uint32_t) (abo_d_tmp >= 2);
        
        //Age
        uint32_t age_d = std::stoi(donor_attributes[2]);
        
        //Sex
        uint32_t sex_d = std::stoi(donor_attributes[3]);
        
        //Size
        uint32_t size_d = std::stoi(donor_attributes[4]);

        //std::cout << "Parsed Donor" << std::endl;

        // Parsing Recipient Attributes
        // HLA
        uint32_t* hla_r = new uint32_t[n_hla];
        counter = 0;
        std::istringstream hla_r_stream(recipient_attributes[0]); 
        while(std::getline(hla_r_stream, tmp, ' ')) {
            if(tmp.empty()) {
                continue;
            }
            hla_r[counter++] = std::stoi(tmp); 
        }

        // HLA_A
        uint32_t* hla_a_r = new uint32_t[n_hla];
        counter = 0;
        std::istringstream hla_a_r_stream(recipient_attributes[1]); 
        while(std::getline(hla_a_r_stream, tmp, ' ')) {
            if(tmp.empty()) {
                continue;
            }
            hla_a_r[counter++] = std::stoi(tmp); 
        }

        //ABO
        uint32_t abo_r_tmp = std::stoi(recipient_attributes[2]);
        uint32_t* abo_r = new uint32_t[2];
        abo_r[0] = (uint32_t) (abo_r_tmp % 2 == 1);
        abo_r[1] = (uint32_t) (abo_r_tmp >= 2);
        
        //Age
        uint32_t age_r = std::stoi(recipient_attributes[3]);
        
        //Sex
        uint32_t sex_r = std::stoi(recipient_attributes[4]);
        
        //Size
        uint32_t size_r = std::stoi(recipient_attributes[5]);

        //std::cout << "Parsed Recipient" << std::endl;
        
        std::shared_ptr<PlainPair> pair = std::make_shared<PlainPair>(n_hla);
        pair->setDonorAttributes(hla_d, abo_d, age_d, sex_d, size_d);
        pair->setRecipientAttributes(hla_r, hla_a_r, abo_r, age_r, sex_r, size_r);
        pairs.push_back(pair);
    }
    
    return pairs;
}


uint32_t evalCrossmatch(uint32_t *hla_d, uint32_t *hla_a_r, uint32_t n_hla) {
    uint32_t *matches = new uint32_t[n_hla];

    for(size_t i = 0; i < n_hla; ++i) {
        matches[i] = hla_d[i] & hla_a_r[i];
    } 

    uint32_t res = 0;
    for(size_t i = 0; i < n_hla; ++i) {
        res = res | matches[i];
    }

    res = ~res;
    res = 1 & res;

    delete[] matches;

    return res;
}


//Update in thesis-draft -> MUX gate at the end is wrong
uint32_t evalHLAmatch(uint32_t *hla_d, uint32_t *hla_r, uint32_t n_hla) {
    uint32_t *mismatches = new uint32_t[n_hla];

    for(size_t i = 0; i < n_hla; ++i) {
        mismatches[i] = hla_d[i] ^ hla_r[i];
    }

    uint32_t sum = 0;
    for(size_t i = 0; i < n_hla; ++i) {
        sum += mismatches[i];
    }

    delete[] mismatches;

    // std::cout << "HLA sum:\t" << sum << std::endl;
    return sum == 0 ? 5 : (sum < 3 ? 3 : (sum < 5 ? 1 : 0));
}


uint32_t evalABOmatch(uint32_t *abo_d, uint32_t *abo_r) {
    uint32_t resL;

    resL = (abo_d[0] ^ abo_r[0]) | (abo_d[1] ^ abo_r[1]);
    resL = ~resL;
    resL = resL & 1;

    uint32_t resR, tmp1, tmp2;
    tmp1 = ~abo_d[0];
    tmp1 = tmp1 & 1;
    tmp2 = ~abo_d[1];
    tmp2 = tmp2 & 1;
    resR = (tmp1 & abo_r[1]) | (tmp2 & abo_r[0]);

    // std::cout << "ABO:\t" << ((bool) (resL | resR) ? 3 : 0) << std::endl;
    return (bool) (resL | resR) ? 3 : 0;
}


uint32_t evalAgematch(uint32_t age_d, uint32_t age_r) {
    bool same = age_d == age_r;
    uint32_t tmp = ~age_d;
    tmp = tmp & 1;
    bool ydor = tmp & age_r;
    
    // std::cout << "Age:\t" << (ydor ? (same ? 2 : 1) : (same ? 2 : 0)) << std::endl;
    return ydor ? (same ? 2 : 1) : (same ? 2 : 0);
}


uint32_t evalSexmatch(uint32_t sex_d, uint32_t sex_r) {
    bool same = sex_d == sex_r;
    uint32_t tmp = ~sex_r;
    tmp = tmp & 1;
    bool fdmr = tmp & sex_d;

    // std::cout << "Sex:\t" << (fdmr ? (same ? 2 : 0) : (same ? 2 : 1)) << std::endl; 
    return fdmr ? (same ? 2 : 0) : (same ? 2 : 1);
}


uint32_t evalSizematch(uint32_t size_d, uint32_t size_r) {
    // std::cout << "Size:\t" << (size_d < size_r ? 0 : 1) << std::endl;
    return size_d < size_r ? 0 : 1;
}


void writeCompGraphToFilePlain(std::vector<std::vector<uint32_t>> &comp_graph, uint32_t n_pairs) {
    std::ofstream outfile("comp_graph_plain.txt");
    outfile << "#CompatibilityGraphPlain" << std::endl;
    for(size_t i = 0; i < n_pairs; ++i) {
        for(size_t j = 0; j < n_pairs; ++j) {
            outfile << comp_graph[i][j] << ",";
        }
        outfile << "\n";
    }
}

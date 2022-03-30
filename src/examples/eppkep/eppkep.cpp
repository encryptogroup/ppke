//Utility libs
#include <ENCRYPTO_utils/crypto/crypto.h>
#include <ENCRYPTO_utils/parse_options.h>
//ABY Party class
#include "../../abycore/aby/abyparty.h"

#include "common/matching_circuit.hpp"
#include "common/compute_cycles_circuit.hpp"
#include "common/find_largest_set_circuit.hpp"

int32_t read_test_options(int32_t* argcp, char*** argvp, e_role* role,
		uint32_t* bitlen, uint32_t* nvals, uint32_t* secparam, std::string* address,
		uint16_t* port, int32_t* test_op, std::string* path, bool* factors,
		uint32_t* n_pairs, uint32_t* n_hla, uint32_t* cycle_length) {

	uint32_t int_role = 0, int_port = 0, bool_factors = 0, int_n_pairs = 0, int_n_hla = 0, 
			int_cycle_length = 0;

	parsing_ctx options[] =
			{ { (void*) &int_role, T_NUM, "r", "Role: 0/1", true, false }, 
					{ (void*) nvals, T_NUM, "nv",
					"Number of parallel operation elements", false, false }, 
					{ (void*) bitlen, T_NUM, "b", "Bit-length, default 32", false,
					false }, 
					{ (void*) secparam, T_NUM, "s", 
					"Symmetric Security Bits, default: 128", false, false }, 
					{ (void*) address, T_STR, "a", "IP-address, default: localhost", 
					false, false }, 
					{ (void*) &int_port, T_NUM, "p", "Port, default: 7766", false,
					false }, 
					{ (void*) test_op, T_NUM, "t",
					"Single test (leave out for all operations), default: off",
					false, false },
					{ (void*) path, T_STR, "d", 
					"Path pointing to the location where the data is stored. Default: '../data/'", 
					false, false},  
					{ (void*) &bool_factors, T_NUM, "f", "Factors to be included 0/1, default: 1",
					false, false},   
					{ (void*) &int_n_pairs, T_NUM, "x", "Number of participating pairs. Default: 10",
					false, false},
					{ (void*) &int_n_hla, T_NUM, "y", "Number of HLA used. Default 50", false, false},
                    { (void*) &int_cycle_length, T_NUM, "z", "Length of cycles. Default: 3", false, false},
					};

	if (!parse_options(argcp, argvp, options,
			sizeof(options) / sizeof(parsing_ctx))) {
		print_usage(*argvp[0], options, sizeof(options) / sizeof(parsing_ctx));
		std::cout << "Exiting" << std::endl;
		exit(0);
	}

	assert(int_role < 2);
	*role = (e_role) int_role;

	if (int_port != 0) {
		assert(int_port < 1 << (sizeof(uint16_t) * 8));
		*port = (uint16_t) int_port;
	}

	*factors = (bool) bool_factors;
	
	if(int_n_pairs != 0) {
		assert(int_n_pairs > 1);
		*n_pairs = int_n_pairs;
	}

	if(int_n_hla != 0) {
		assert(int_n_hla > 0);
		*n_hla = int_n_pairs;
	}

    if(int_cycle_length != 0) {
		assert(int_cycle_length > 1);
		*cycle_length = int_cycle_length;
	}

	//delete options;

	return 1;
}


int main(int argc, char** argv) {
	
	e_role role;
	uint32_t bitlen = 32, nvals = 8, secparam = 128, nthreads = 1, n_pairs = 10, n_hla = 50, cycle_length = 3;
	uint16_t port = 7766;
	std::string address = "127.0.0.1";
	std::string path = "../data/input/data_1000/";
	int32_t test_op = -1;
	e_mt_gen_alg mt_alg = MT_OT;
	bool factors = true;

	read_test_options(&argc, &argv, &role, &bitlen, &nvals, &secparam, &address,
			&port, &test_op, &path, &factors, &n_pairs, &n_hla, &cycle_length);

	seclvl seclvl = get_sec_lvl(secparam);

	std::cout << "Role " << role << std::endl;

	std::cout << "Part 1 Matching...\n" << std::endl;

	matching_circuit(role, address, port, seclvl, bitlen,
			nthreads, mt_alg, path, factors, n_pairs, n_hla);

	std::cout << "\n\n##########" << std::endl;
	std::cout << "##########\n" << std::endl;

	std::cout << "\nPart 2 Compute number of cycles...\n" << std::endl;

    auto n_cycles = compute_cycles_circuit(role, address, port, seclvl, bitlen,
			                nthreads, mt_alg, n_pairs, cycle_length);

	std::cout << "\n\n##########" << std::endl;
	std::cout << "##########\n" << std::endl;

	std::cout << "\nPart 3 find solution...\n" << std::endl;

    find_largest_set_circuit(role, address, port, seclvl, bitlen,
			nthreads, mt_alg, n_pairs, cycle_length, n_cycles);

	return 0;
}


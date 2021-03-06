//Utility libs
#include <ENCRYPTO_utils/crypto/crypto.h>
#include <ENCRYPTO_utils/parse_options.h>
//ABY Party class
#include "../../abycore/aby/abyparty.h"

#include "common/find_largest_set_circuit.hpp"

int32_t read_test_options(int32_t* argcp, char*** argvp, e_role* role,
		uint32_t* bitlen, uint32_t* nvals, uint32_t* secparam, std::string* address,
		uint16_t* port, int32_t* test_op, uint32_t* n_pairs, uint32_t* cycle_length) {

	uint32_t int_role = 0, int_port = 0, int_n_pairs = 0, int_cycle_length = 0;

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
					{ (void*) &int_n_pairs, T_NUM, "x", "Number of participating pairs, default: 10",
					false, false},
					{ (void*) &int_cycle_length, T_NUM, "y", "Length of cycles, default: 3", false, false}
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

	if(int_n_pairs != 0) {
		assert(int_n_pairs > 1);
		*n_pairs = int_n_pairs;
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
	uint32_t bitlen = 32, nvals = 8, secparam = 128, nthreads = 1, n_pairs = 10, cycle_length = 2;
	uint16_t port = 7766;
	std::string address = "127.0.0.1";
	int32_t test_op = -1;
	e_mt_gen_alg mt_alg = MT_OT;

	read_test_options(&argc, &argv, &role, &bitlen, &nvals, &secparam, &address,
			&port, &test_op, &n_pairs, &cycle_length);
			
	seclvl seclvl = get_sec_lvl(secparam);

	uint32_t n_cycles = 18;    	// 27 cycles for 5 vertice cycles length 3
                                // 10 = 5 vertices cycle length 2

                                // data_100 files
                                // xxx = 10 vertices cycle length 2
                                // yyy = 10 vertices cycle length 3
                                // 128 = 20 vertices cycle length 2
                                // 387 = 20 vertices cycle length 3
                                // 1107 = 50 vertices cycle length 3
                                // 2307 = 100 vertices cycle length 3

	std::cout << "Role " << role << std::endl;

	find_largest_set_circuit(role, address, port, seclvl, bitlen,
			nthreads, mt_alg, n_pairs, cycle_length, n_cycles);

	return 0;
}


#ifndef __PAIR__HPP__
#define __PAIR__HPP__

#include "patient.hpp"

#include <memory>

class share;
using share_p = std::shared_ptr<share>;

class Pair {

public:
    //Pair(share_p id, std::shared_ptr<Patient> donor, std::shared_ptr<Recipient> recipient);
    Pair(int id, std::shared_ptr<Patient> donor, std::shared_ptr<Recipient> recipient);
    ~Pair();

	//share_p getId() const;
    int getId() const;
    std::shared_ptr<Patient>& getDonor();
    std::shared_ptr<Recipient>& getRecipient();

	
private:
	//share_p id;
    int id;
    std::shared_ptr<Patient> donor;
    std::shared_ptr<Recipient> recipient;
};

#endif
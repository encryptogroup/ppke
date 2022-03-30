#include "pair.hpp"

// #include "CircuitWrapper.h"

// Pair::Pair(share_p id, std::shared_ptr<Patient> donor, std::shared_ptr<Recipient> recipient) :
// 	id(id), donor(donor), recipient(recipient) {

// }

Pair::Pair(int id, std::shared_ptr<Patient> donor, std::shared_ptr<Recipient> recipient) :
	id(id), donor(donor), recipient(recipient) {

}

Pair::~Pair() {

}

// share_p Pair::getId() const {
//     return this->id;
// }

int Pair::getId() const {
    return this->id;
}

std::shared_ptr<Patient>& Pair::getDonor() {
    return this->donor;
}

std::shared_ptr<Recipient>& Pair::getRecipient() {
    return this->recipient;
}
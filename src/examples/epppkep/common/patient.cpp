#include "patient.hpp"

#include "CircuitWrapper.h"

Patient::Patient(share_p hla, std::vector<share_p> bloodGroup, share_p age, share_p sex, share_p size) :
	hla(hla), bloodGroup(bloodGroup), age(age), sex(sex), size(size) {
}

Patient::~Patient() {    
}

share_p& Patient::getHLA() {
    return this->hla;
}


std::vector<share_p>& Patient::getBloodGroup() {
    return this->bloodGroup;
}

share_p& Patient::getAge() {
    return this->age;
}

share_p& Patient::getSex() {
    return this->sex;
}

share_p& Patient::getSize() {
    return this->size;
}

Recipient::Recipient(share_p hla, share_p hla_a, std::vector<share_p> bloodGroup, share_p age, share_p sex, share_p size) :
    Patient(hla, bloodGroup, age, sex, size), hla_a(hla_a) {
        
}

Recipient::~Recipient() {

}

share_p& Recipient::getHLA_A() {
	return this->hla_a;
}
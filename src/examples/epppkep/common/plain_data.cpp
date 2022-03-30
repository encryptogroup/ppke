#include "plain_data.hpp"

PlainPair::PlainPair(uint32_t n_hla) :
    n_hla(n_hla) {
}

PlainPair::~PlainPair() {
    delete[] hla_d;
    delete[] abo_d;
    delete[] hla_r;
    delete[] hla_a_r;
    delete[] abo_r;
}

uint32_t PlainPair::getN_HLA() const {
    return this->n_hla;
}

void PlainPair::setDonorAttributes(uint32_t* hla, uint32_t* abo, uint32_t age, 
                                    uint32_t sex, uint32_t size) {
    this->hla_d = hla;
    this->abo_d = abo;
    this->age_d = age;
    this->sex_d = sex;
    this->size_d = size;
}

void PlainPair::setRecipientAttributes(uint32_t* hla, uint32_t* hla_a, uint32_t* abo, 
                                        uint32_t age, uint32_t sex, uint32_t size) {
    this->hla_r = hla;
    this->hla_a_r = hla_a;
    this->abo_r = abo;
    this->age_r = age;
    this->sex_r = sex;
    this->size_r = size;
}

uint32_t* PlainPair::getHLA_D() {
    return this->hla_d;
}

uint32_t* PlainPair::getABO_D() {
    return this->abo_d;
}

uint32_t PlainPair::getAge_D() const {
    return this->age_d;
}

uint32_t PlainPair::getSex_D() const {
    return this->sex_d;
}

uint32_t PlainPair::getSize_D() const {
    return this->size_d;
}

uint32_t* PlainPair::getHLA_R() {
    return this->hla_r;
}

uint32_t* PlainPair::getHLA_A_R() {
    return this->hla_a_r;
}

uint32_t* PlainPair::getABO_R() {
    return this->abo_r;
}

uint32_t PlainPair::getAge_R() const {
    return this->age_r;
}

uint32_t PlainPair::getSex_R() const {
    return this->sex_r;
}

uint32_t PlainPair::getSize_R() const {
    return this->size_r;
}
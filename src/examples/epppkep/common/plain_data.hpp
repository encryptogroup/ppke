#ifndef __PLAIN_DATA__HPP__
#define __PLAIN_DATA__HPP__

#include <memory>

class PlainPair {

public:
    PlainPair(uint32_t n_hla);

    ~PlainPair();

    uint32_t getN_HLA() const;

    void setDonorAttributes(uint32_t* hla, uint32_t* abo, uint32_t age, 
                            uint32_t sex, uint32_t size);
    void setRecipientAttributes(uint32_t* hla, uint32_t* hla_a, uint32_t* abo, 
                                uint32_t age, uint32_t sex, uint32_t size);

    uint32_t* getHLA_D();
    uint32_t* getABO_D();
    uint32_t getAge_D() const;
    uint32_t getSex_D() const;
    uint32_t getSize_D() const;

    uint32_t* getHLA_R();
    uint32_t* getHLA_A_R();
    uint32_t* getABO_R();
    uint32_t getAge_R() const;
    uint32_t getSex_R() const;
    uint32_t getSize_R() const;

private:
    uint32_t n_hla;

    uint32_t* hla_d;
    uint32_t* abo_d;
    uint32_t age_d;
    uint32_t sex_d;
    uint32_t size_d;

    uint32_t* hla_r;
    uint32_t* hla_a_r;
    uint32_t* abo_r;
    uint32_t age_r;
    uint32_t sex_r;
    uint32_t size_r;

};

#endif
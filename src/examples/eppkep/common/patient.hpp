#ifndef __PATIENT__HPP__
#define __PATIENT__HPP__

#include <memory>
#include <vector>

class share;
using share_p = std::shared_ptr<share>;

class Patient {
public:
	Patient(share_p hla, std::vector<share_p> bloodGroup, share_p age, share_p sex, share_p size);

    virtual ~Patient();

    share_p& getHLA(); // Dont pass by reference because we do not want the object to be changed once it is set, and passing by const reference causes significant casting overhead for using the gates
    std::vector<share_p>& getBloodGroup();
    share_p& getAge();
    share_p& getSex();
    share_p& getSize();

private:
    share_p hla; // Contains all HLA in the following Order -> HLA-A, -B, -DR, and -DQ
    std::vector<share_p> bloodGroup; // 2 Bit variable
    share_p age;  //True = junior (<50 years old), False = senior(>50 years old)
    share_p sex; //True = Female, False = Male
    share_p size; // rounded to the next integer
};

class Recipient : public Patient {

public:
    Recipient(share_p hla, share_p hla_a, std::vector<share_p> bloodGroup, share_p age, share_p sex, share_p size);

    virtual ~Recipient();

    share_p& getHLA_A();

private:
    share_p hla_a; // HLA Antibodies -> HLA-A, -B, -DR, and -DQ
};

#endif
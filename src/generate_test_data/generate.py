import generate_data as gd

def write_pair_to_file(pair, nb: int) -> None:
    
    name = "../../data_100/pair" + str(nb) + ".txt"
    f = open(name, "w")
    f.write("# Number of HLA\n")
    f.write(str(50) +" \n\n")

    donor = pair[0]
    f.write("### Donor attributes\n")
    f.write("# HLA\n")
    for i in donor[0]: 
        f.write(str(i) + " ")
    f.write("\n")
    f.write("# ABO\n")
    f.write(str(donor[1])  + "\n")
    f.write("# Age\n")
    f.write(str(donor[2])  + "\n")
    f.write("# Sex\n")
    f.write(str(donor[3])  + "\n")
    f.write("# Size\n")
    f.write(str(donor[4])  + "\n\n")

    f.write("#"* 10 + "\n\n")
    recipient = pair[1]
    f.write("### Recipient attributes\n")
    f.write("# HLA\n")
    for i in recipient[0]: 
        f.write(str(i) + " ")
    f.write("\n")
    f.write("# HLA_A\n")
    for i in recipient[1]: 
        f.write(str(i) + " ")
    f.write("\n")
    f.write("# ABO\n")
    f.write(str(recipient[2])  + "\n")
    f.write("# Age\n")
    f.write(str(recipient[3])  + "\n")
    f.write("# Sex\n")
    f.write(str(recipient[4])  + "\n")
    f.write("# Size\n")
    f.write(str(recipient[5])  + "\n")


for i in range(100):
    pair = gd.DataGenerator(i * 101).generate_pair()
    write_pair_to_file(pair, i)



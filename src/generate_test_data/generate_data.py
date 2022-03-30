import random
import numpy as np

class DataGenerator(object):

    def __init__(self, start) -> None:
        self.seed = start

    def generate_hla(self) -> np.array:
        hla = np.full(shape = 50,
            fill_value = 0,
            dtype=int)

        # Generate 6 unique indices to set them to true
        entries = random.sample(range(0, 50), 6)
        for x in entries:
            hla[x] = 1
        
        return hla
    
    def generate_hla_a(sel, hla: np.array) ->np.array:
        hla_a = np.full(shape = 50,
            fill_value = 0,
            dtype=int)

        # select a random index from the positive entries in the given hla array
        entries = random.sample(range(0, 50), 2)
        for x in entries:
            hla_a[x] = 1
        
        return hla_a

    def generate_abo(self) -> int:
        return random.randint(0,3)

    def generate_age(self) -> bool:
        return random.randint(0, 1)

    def generate_sex(self) -> bool:
        return random.randint(0,1)

    def generate_size(self) -> int:
        return random.randint(60, 90)

    def generate_pair(self):
        self.seed += (int) (random.random() * self.seed)
        random.seed(self.seed)
        donor = (self.generate_hla(), self.generate_abo(), self.generate_age(), self.generate_sex(), self.generate_size())
        recipient = ((self.generate_hla(), self.generate_hla_a(donor[0]),self.generate_abo(), self.generate_age(), self.generate_sex(), self.generate_size()))

        return (donor, recipient)
    
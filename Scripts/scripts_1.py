import Bio
from Bio.Seq import Seq


#sequence = "ACAATGAGGTCACTATGTTCGAGCTCTTCAAACCGGCTGCGCATACGCAGCGGCTGCCATCCGATAAGGTGGACAGCGTCTATTCACGCCTTCGTTGGCAACTTTTCATCGGTATTTTTGTTGGCTATGCAGGCTACTATTTGGTTCGTAAGAACTTTAGCTTGGCAATGCCTTACCTGATTGAACAAGGCTTTAGTCGTGGCGATCTGGGTG"
def find_start_codons(sequence):
    """
    Find all start codons in a given DNA sequence and its complementary strand and return their positions.

    Args:
        sequence (str): The DNA sequence to search for start codons.

    Returns:
        dict: A dictionary containing the positions of all start codons found in the sequence.

    """

    #if sequence.endswith(".txt"):
    with open(sequence, "r") as file:
        sequence = file.read()

    start_codon = "ATG"
    sequence = Seq(sequence)
    reverse_complement = sequence.reverse_complement()
    forward_codon_dict = {}
    reverse_codon_dict = {}

    for i in range(len(sequence) - 2):
        if sequence[i:i+3] == start_codon:
            forward_codon_dict[i] = sequence[i:i+3]

    for i in range(len(reverse_complement) - 2):
        if reverse_complement[i:i+3] == start_codon:
            reverse_codon_dict[len(sequence) - i - 3] = reverse_complement[i:i+3]
    file.close()
    return forward_codon_dict, reverse_codon_dict

#find_start_codons(sequence)

def main():
    forward_codon_dict, reverse_codon_dict= find_start_codons("Bioinformatics_Portfolio/Data/Vibrio_cholera_genome.txt")
    print("forward_codon_dict",forward_codon_dict)
    print("reverse_codon_dict",reverse_codon_dict)
    

if __name__ == "__main__":
    main()

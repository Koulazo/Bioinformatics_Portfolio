import Bio
from Bio.Seq import Seq

def kmer_detector(Text,k):
    """
    Find all kmers in a given DNA sequence and its complementary strand and return their positions.

    Args:
        sequence (str): The DNA sequence to search for start codons.

    Returns:
        dict: A dictionary containing the positions of all start codons found in the sequence.

    """

    #if sequence.endswith(".txt"):
    # with open(Text, "r") as file:
    #     Text = file.read()
    highest_kmer = 0
    kmer_dict = {}
    for i in range(len(Text) - int(k) + 1):
        kmer = Text[i:i+int(k)]
        if kmer not in kmer_dict:
            kmer_dict[kmer] = 1
        else:
            kmer_dict[kmer] += 1
    for i in kmer_dict:
        # if kmer_dict[i] == 1:
        #     pass
        if kmer_dict[i] > highest_kmer:
            highest_kmer = kmer_dict[i]
        elif kmer_dict[i] < highest_kmer:
            pass
    print("kmer_dict",kmer_dict)
    print("highest_kmer",highest_kmer)
    return kmer_dict


def reverse_complement(sequence):
    """
    Find all start codons in a given DNA sequence and its complementary strand and return their positions.

    Args:
        sequence (str): The DNA sequence to search for start codons.

    Returns:
        dict: A dictionary containing the positions of all start codons found in the sequence.

    """
    sequence = Seq(sequence)
    reverse_complement = sequence.reverse_complement()
    print("reverse_complement:",reverse_complement)
    return reverse_complement


def find_occurrences(sequence, pattern):
    """
    Find all occurrences, forwards and backwards, of a pattern within a nucleotide sequence.

    Args:
        sequence (str): The nucleotide sequence to search for occurrences.
        pattern (str): The pattern to search for.

    Returns:
        list: A list of all occurrences found in the sequence.

    """
    occurrences = []
    sequence_length = len(sequence)
    pattern_length = len(pattern)

    for i in range(sequence_length - pattern_length + 1):
        if sequence[i:i+pattern_length] == pattern:
            occurrences.append(i)
        if sequence[i:i+pattern_length] == pattern[::-1]:
            occurrences.append(i)
    print("occurrences",occurrences)
    return occurrences


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
    #file.close()
    return forward_codon_dict, reverse_codon_dict


def main():
    Text = "TCAAGTCATGGGACAATCTCAAGTCATGTCAAGTCATGGTAGTGTAGCTTGATTCACATTCTCCCTGGACAATCTCAAGTCATGTTGATTCAGTAGTGTAGCTTGATTCATTGATTCAGGACAATCTTGATTCATTGATTCATCAAGTCATGTTGATTCACATTCTCCCTTTGATTCAGGACAATCTTGATTCAGGACAATCGGACAATCGGACAATCCATTCTCCCTGTAGTGTAGCTTGATTCATCAAGTCATGGGACAATCCATTCTCCCTTTGATTCAGTAGTGTAGCCATTCTCCCTTTGATTCAGGACAATCCATTCTCCCTTCAAGTCATGTTGATTCAGGACAATCTTGATTCAGTAGTGTAGCCATTCTCCCTCATTCTCCCTCATTCTCCCTTTGATTCAGGACAATCCATTCTCCCTCATTCTCCCTTCAAGTCATGTTGATTCAGGACAATCTTGATTCATTGATTCAGGACAATCCATTCTCCCTTCAAGTCATGTCAAGTCATGTCAAGTCATGCATTCTCCCTCATTCTCCCTGGACAATCGGACAATCCATTCTCCCTTTGATTCATTGATTCACATTCTCCCTCATTCTCCCTCATTCTCCCTGGACAATCTTGATTCACATTCTCCCTTTGATTCATCAAGTCATGGTAGTGTAGCGGACAATCCATTCTCCCTCATTCTCCCTCATTCTCCCTGGACAATCGGACAATCGGACAATCGGACAATCGGACAATCTCAAGTCATGCATTCTCCCTGTAGTGTAGCTTGATTCATTGATTCAGTAGTGTAGCTTGATTCAGGACAATCGGACAATCGGACAATCTCAAGTCATGTCAAGTCATGGGACAATCGTAGTGTAGCTCAAGTCATGGTAGTGTAGCGGACAATCGGACAATCTTGATTCATTGATTCATCAAGTCATGCATTCTCCCTCATTCTCCCTGTAGTGTAGC"
    pattern = "ATAT"
    k= "12"
    sequence = "GATATATGCATATACTT"
    forward_codon_dict, reverse_codon_dict= find_start_codons("/home/namorak/Bioinformatics_Portfolio/Data/Vibrio_cholera_genome.txt")
    print("forward_codon_dict",forward_codon_dict)
    #kmer_detector(Text,k)
    #reverse_complement(sequence)
    find_occurrences(sequence, pattern)

if __name__ == "__main__":
    main()
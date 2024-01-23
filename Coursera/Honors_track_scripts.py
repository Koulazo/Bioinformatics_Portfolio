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




def main():
    Text = "TCAAGTCATGGGACAATCTCAAGTCATGTCAAGTCATGGTAGTGTAGCTTGATTCACATTCTCCCTGGACAATCTCAAGTCATGTTGATTCAGTAGTGTAGCTTGATTCATTGATTCAGGACAATCTTGATTCATTGATTCATCAAGTCATGTTGATTCACATTCTCCCTTTGATTCAGGACAATCTTGATTCAGGACAATCGGACAATCGGACAATCCATTCTCCCTGTAGTGTAGCTTGATTCATCAAGTCATGGGACAATCCATTCTCCCTTTGATTCAGTAGTGTAGCCATTCTCCCTTTGATTCAGGACAATCCATTCTCCCTTCAAGTCATGTTGATTCAGGACAATCTTGATTCAGTAGTGTAGCCATTCTCCCTCATTCTCCCTCATTCTCCCTTTGATTCAGGACAATCCATTCTCCCTCATTCTCCCTTCAAGTCATGTTGATTCAGGACAATCTTGATTCATTGATTCAGGACAATCCATTCTCCCTTCAAGTCATGTCAAGTCATGTCAAGTCATGCATTCTCCCTCATTCTCCCTGGACAATCGGACAATCCATTCTCCCTTTGATTCATTGATTCACATTCTCCCTCATTCTCCCTCATTCTCCCTGGACAATCTTGATTCACATTCTCCCTTTGATTCATCAAGTCATGGTAGTGTAGCGGACAATCCATTCTCCCTCATTCTCCCTCATTCTCCCTGGACAATCGGACAATCGGACAATCGGACAATCGGACAATCTCAAGTCATGCATTCTCCCTGTAGTGTAGCTTGATTCATTGATTCAGTAGTGTAGCTTGATTCAGGACAATCGGACAATCGGACAATCTCAAGTCATGTCAAGTCATGGGACAATCGTAGTGTAGCTCAAGTCATGGTAGTGTAGCGGACAATCGGACAATCTTGATTCATTGATTCATCAAGTCATGCATTCTCCCTCATTCTCCCTGTAGTGTAGC"
    Pattern = "ATG"
    k= "12"
    sequence = "GGACAGTGCAGTTAGCCGAATGGTTTACATCCTGTGAAAGGTGCCACTAGCGGCTGATCAGCGGACTCCCTAAAAAACCATCTTTGCATTTGACGCGCCATGGCGGACGCTTCGGCTTATAGTGGATTCAAATTATCTTGGCTGCCCAACCGGACCATACCACCCCTTCGCCAGTGACAAAACTACCCCGGGCGCGTTTCGTTCAGCCTGGGGGTCAAAACAGACGACCGACAATATTCGAATACAAGTAGTTCGCCGAACCGGCTAGGCGTACAGCAGTCTATGCTGCACTCTCGGATTGCCTAAAGGTTAGGGAATCTGACAGTACTGAACGTCTACCTCCTGCCCCTTCCCGAATGCGTGTGTATGTTCTTGTGATAGCTGCACTGACGTTACTTAGCACCGTATAAGACAGTTTGTTCGAAGTATGCAGTTACCGCTGTGATTTTCCCTCTCAAGTAAGGCTGCACGCCAATGCCAAGAGTTACTGGCATGCACCTACACCGGATAGAGTACACACAGAAACTTTCTTAGCAGACCAATTTATCAGATTCATTCGTATGGGCCCACGACAGCAGAAAATTTGGCCCTAAGTACGCTTTCCCCGCCAGGGGCTCGATGGCGGCCGCGTTACCGGTGACATTTGGTACTCAGAACTCCACCGCATCGTGCGCGTTGGATAGACTAAAAATAGCATTGATCCAGTCTCTGGGTCGCGCAGCCCGTTGAGACAAACGACGTATGGAGTACTGTCGGCGAAGATAAGTCTTTTTTACGCGGGAGTCGTCAGCACTGACGTACAACAAACATAAGCTATTCACCTGACTGGAGCCGTCTCGCTGTGGCGGGAGGTGCGGCAACCGAACGGGACGTGTGACAGCATGTCCTAATCTAACTACCATATTATACGAGCTCGGCTCAACGGTCAGGAAACATAGATAGGGATGGCAGATGCTATTCATCCCATCCAGTCTCATATATCGCCTCATGGCCATATGATTCCTTCGCTCACCCTTCTATTCCCCAAGATGCGTCAGTTTGAGGTTAGGGGAGCGTTTAAGTCCATGGCGGATACCTCGTCCTGCCTGTGTTATCTTCCTGACCTAGCGGGGCTCTCGCGATATTCCCCTCGCGGTCACAGACTTTGCGAATTAGGGGGGCCGGCGTTGCTTAAACACACAGCGTTTCCCGAACCAGGACATCACGCCTAGGTCCGTCATACGACTTCGGTGTACGAACCCTAGACCTAATAACCGCTCACGATCTCGTGATCAGATCATTCGCGCTAAGACCTCTTTCCCGAATTTGCATGAGTGTGTTACTCCGGTGCGGATCCAGCCGCCACCCAATCGGAACATAGCTGATCGCTACGTTAGTCGGAGCGTGGGGATTAAACTCTTTAGGGTTCAATGGGAGTAGTCTCAAACTAGAGCGGGTTGTTCGTAGATAATCTTGCCACATCATGCGGTGTTGTCTGGCTCAACTCGATTTCTAGCACCATCAAAACGAAACAATGATCTTTAAGATGACTAACTCTTCCATGATCGTATTGCTAAGCAACTTCGGTTGCGCGCGTGTGATGCGCAAGATCGGCCGTCCTTCACCGATCTTTGATACGAGGCATGCGTAGCCTCTTGAAGCGCTGGGTCAGGCTGCTCGGTCCTTCCCGAAGTCTGCAAGAACCATACAAGGTTTGTCTAATCAGTTATACTCCCTATCCAGTGTCGAATGGGAGATCATTTGCGATTGTGTCACTCGCATGGTGCCAGTATATCACAGGACAGCCCTATTTATACCTCGTAGCTGTTAAGGGACGACTCAAACTCTGCACATCACTTTAGGGATGTTGAAATATGTGACGTTCTCACCGCTGCAGTATGCAAACCCCGACATTAGGCCGTCCATAGTAATGGCGGGATCGGAAGGGCCCAGAGTCACGGTTGTCGTGAGTAGCCTTTTAACCGGACGCTGTGGAAGTAAAAATTTCGCACACTTATAGCTCAAGATCCCAAGAATGGGGACTGAGGGTACCTTTTCGACCAGGGCAACAAGTCCGGCGTCTCTCCAGGACAAACGATTCTTATGGCCTGCCACCATTCGATCCCGGTTATCGTGTCCAGCGTTCCCTCAGATCTCGCATGTGCGGTCTCTTCTCAAGGTTCTGGGGCCAAAGAGTCATTGCGTACCCACGAGGTGCGACACAGTAGAGAATGTAGTAAAGCATATTACTCGGGGCGCTTTGGTTAAAACACCAGCGGGGGTCCAGATAAACTACGGAGCAACACCGGTGAACGAATCACCGCGTCCCCGTATGTGTCGGTGGCGCGAGACTAATCATGCTTGATACATGGAAGCTGAGTTGACTTACAATCCTATCAGGTGGTGAACGGTCTAAGCAATTGATTCCAGTGGCTTGTTATAAGTCTGACGCTTACCTGAAGCTGATCAACAGGTCTCACACACCTATTCGGCTGTTAGAGTAAGAAGACAATGGTGCCTTCAGTTTCCCTATCAATACCGACCCCCCCCTTTGTCAATACCGCGGAGAATAGACGGCGTGAGTCAGGCAAGGCCTGTATCCTGGTGCGCAGTGCCTCTTGACTGGCACTGGCGAATCCTTAGTAAGCAGTCCCCAGATAAGGAAATCTTTCATTTATGTAACAGTCAGTTCGCCGTGGATAGAGACGGCAGTTCACGCGCTGCATTAACCCATGCATTAGAAAACTCGTCGATTGCCGGATGCCGCGGAGTGAATGACAGAATCACGGAATAATAATCGATCATGGCATTAGGTCGTGTAGGGCACAGCCGACGTCTGATCTTGTGGGCTGCACTTCGTTGCGTAAGTGGATTAGTTAGTAATATTATCAGTTGGAGTATCGGTAATGGCTATTTATCTGATCAATAGTGATCGTTACTATAAACCCACGTCGCGTCTCGCTATTGTTATTCGTTCAGAGCGCGGATATTGTGACCTCTTCACGGTCCCCGTCAGTCATCATACACTCTAACAGTGGTTTTTGGTCTAAGGCTGGGTCATTGGCCCGGTGGGTTCGTGGAGGAGTCAGTCCATCATCGCGGTGCAAAATGAAACTCCATGTGGCCTCGGAAGCCCTCTCCACGCCCAGTCACAGATTAGCCTAGGGAGCATACCGACCAGGACCTACTGATACGCCACGTCACTATACCGTAGATTTAGCATGAGAGCAGATTCCCACGCAGATGTTTTAAACTTTAATTTAGGTACTTACACTCTTATATGTTAATTTCCTGAATGTATGCGTAATACCCAGTAGGATAGGCCTCTCGGTGGCTTGTTCAACTGTATCTAAGCTTACGACCAGACTACAAGGTGTCAGGATGTACAGCCTCCCGTAAGGAACCGCCCAACACCGCATCGCGCCTCGGTGTGCTGATCTCCGTCTTCCCAACAATCGGTCGAGGATGACTGATGCGGTGGTATGGTGTTGATGTAAATCTAGAGACGAACTACAATTGCCCCAATCCCGGAGAGCGGGTAAGCGTAAAGAGCGGTTTTAGATTTCAAAAATTAAATCCGGATACGCGATCACATATTCGTTATCCTGATTAGGTGTCCGTGAAGTCATAAAAATCGGTACTTGTTCATAGCTGTCATTGATCCCTTACATCATTCCCTTAGCCTAAACCTAATCACGGTACGGGTCGTGCAGTGTACGTTCAATTGATCTGGCTGATTTTTGTAAATCACCAGCGTAATGTTGAGTCGGATGCTCCGGACCTGACGCGGTTTACACATTCACTTGAAGAGCCGTGAGGGTATTGCTCTAACAAAGGTTCCGTTTTTGCAGGGGCCGATATAATGGCAAACCTGGCGGTAGAAGACAGTGAAGTTGTACTTCAATCTCTTGGATGTAATCAACCTGCCAGTCGCATAACATCCGAACGAAGAAGGTGTCCGGAGACAATCGTGCAAAAGATAGGTGGCTCGTTCACTTAGCCATGGTAGTATGCGGGGAACGCTAACTTCCTACCTGAGTAACTACTCACGGATTTTCTAAATAAATAACGCCACTGAGTTGTACTAGACTCCTGCAGGAATCGGATTAAGGGCTCCTCCGTTAATAGTTTCTTTACTTTTAAGAAATAATTTGAATCTTCGGTAACCGGACATGCGAAAAGCGATCAATGAACCGCGCTTGCTGGATCTATTCGAGAGGCTGCTGCTCCCAAATCGAAACCGCGGTTACCAGGATAATTGAGCGTATTATGTCTCTTAATTTTCTTCCGGGAAATCACAAGTTTATTTACCGACGCTTCGATTCACTGTTATAATAGGATACTAAACAGCACCAGATTACGGGGTGCCGACCTTAAGCCGTCAAGCTATTCCGGGGTAAGACCATCAAGCCGAGGTAAGCACGGGGTGTGCCCCTAGACCAGAACTCTGCCAACTGGGCAGGCGCGATTCGTTTGTGTGATCTTGAAAACGTTTCTCTGGGAAATATAACACTGGGGAGAAGACATATCAATCGATGATTCACACTAAGAGTATATGGTGAAGCGAGTTTTGTGGGAGGCAATCAACATAAGCATGTTCGTGCGTTTAAAGGTAGTTCACTAAAACAGATCAAGTCCGGCTGGTCCAGATGGAAGCTTGGAGCCATCCTGTCGTGTTCAACTAAGAGTCCTTGGTGTACAGGAAACCCACTCTAATTCCAGCATTCCAGGCATAACATGGCGTGAACCGGCTTACTCATTTGATATGGTCACCTTTTCTCTCAATCTCCATTCGTCTTACGCACTTGGCTCCCCCATCACCGCTACGGGAGCTCCTTTGCGAAGACTTCGCCTGTCGCTTGGCATTCGTAGCTTGCAACTAGCAACTTGAGACGGCTAACACTACAAATATCTGCCCTGCCGGCCTCCCCCCCCTTTGTTCGTGGAATTTACTAAAACTTATCGAGTGGCTTGATGAGTGAGGTTCGTGTGAGCGCCTCGGTGCTGAGTTCTCTGCTAGCTACCTGGCTGCCCTCTCTATTAGTGCCTAATATCCCCAAGTGCTTTTAAGAGACGTCGCGATACTGAGACGTAAGCAGCTTGTGTTCCGTAGCAGGCGCGGAGCTCCCACGGTTAAGATAGGCTTTGGCATGGTCATGCAGCGGATTGACTTGCGCCAGTGAAGAGTACGATCTCTGCTGATATTCCGTCACTGAATTTAAGGTCGGAGGTGCAAGTGCCGGTATCATAAGATATAGTGTGTCGTGAAATGGTTGTTCTACAATGTCACGCGCTCTAAATAGCTCGATATTTAAGTACATGGGACAGTACAGGGCTAGCCTTAGGTGTCCCGCCGTAAACAGCGGCTCTAGACAAGGTTCCGGTCCACTTGATAGGTGCGTAGATCTCCCCAAGGCAATTACCGCATGCGTGCGAGGTTCCCCCAACTGAAGTATGAGCTCGCTACAATATTGGCAACGATCCCACTGTACCCAGCGTCCACCATAGCATAGGTTTATCACATTTTAAATACCAGGCCAAGCGTACAAGTCACTCTTTATTAGTGGCGGCATTGAAATGAATGCATCTCCGCTTGGTCAGGTCATGGCACGAGATTAAGGGTTAGAAGGTTGGTAGCAGCCGTAGTATGCCACGCGCGGCACGCATAGTTCTTTAGTGCCACTGATATCTCGCCACCTGATCAAACGTTGGATCCAGTTACGCAGCTCCAGGTGGTAGTGGCGACGATCCCAGTACCAAACTAAGCGTCTTTCGGTATTCGCAAATTGGGCAAGACTCTCGTAGTGCGCTTGATACCCCTTGCTCAGGCGAGCCACAGTCCGGTGTACTAAAAAGTTATGACGCACCTGTCTGCATATTAATGGCTCCCCTTAGACGCTCCCCTAGGCCATGAGATTTCGCACCGTGCTACCCTATTTAAGCAGGCCGTAAGGTGCTACGGCCACCCGCAAGAAAGGCCACTCTACCGGGCGCCTTCCATCTGCCGTTTCAAATAGACACCTCCAACAAGACACAGCGTCTATGTATTACGGGACTGCTGATCAATGATAAACCGCATGACTTGTGCTTGGTGTTGTAGAGCCAGGAGCCGCCACAATTCTGCACCTAAAGCATTCTTGATGAAATGCTATGGTCGGCACACTATCTGCCGACTACCTCACAAGACTCCAGCTGAAGGGAATCAACCGTACCACATTTCTGTTGAGATGAGTAGTCGCTTTGGTCATGCCGGAACCACTGGTAGAGGGGTCTGCGACAATACACCCCGTCCCTGCAAACTGGATAAGTCGCCTCGAATTCGGCCATGGACGGAAGGGAGTGTGCTGTTGGCCTAGTCCCGCGAATTGACAAGGTAATGAGTGATCCTGCTCTCACTTCCATTGCTAACGGGATCAGCGATCTATAATCACGCGTGCGGCTACCAGACATTGCGCGCAACAAACCGTAAAAGGCTAATAGAGGCCGATCACCAGCTGGCAGGATCATTCTCGAGGACTGTCTAGACGCTGCGTGTCGCACAGCGGCGTGTGAGGGACAGACCCTGGTAGTCGCTTGCTGCCTGGGAACGATATCGTGAAGGTAGTCGAACTGTATATGAATCTCGACCCGTACAAAATGGAACATCGAGGCTTTGTACTACGTGTAGGGGAGGCCTGAGCGCTGCGCCGGTCTCGCTAGGTGCAATGGACACAGATAAGGTCTAAAGCTTACACCCTTTGGGATTCTTCATCAAAGCTCTCGAAATAGCACACAGTATAGAGCGTTGAAGTGCCTTTACCCGACTCCGAGAGCACTACAGATCATTCACGCATATGAGAGTGGACCTATCAGTTAGGTTTACCATGGTCTCGGTATCATGATAGTCGAGAGGTGACACGGCCGGGGGGTTCACAACGGCGAGGGAACGTAGTCATTGCCGTTTATGGCCGTGTTCTAGCCTTCACGGGTGGACGGACGAACTGCAGCAGTGCGTTTGTTTTAAGCAAGTAGCCGCCATAGGGGTCCGATGATCAATTGTATGAGACCCTCGACCTTTTAATTCCTGAGTCACTCCGAAGTACTCTTAGACTCATTTGAGTATCACTTGCAGGTCCCGGTACCAACGATGACCTAAGGAGTTCCTCTCCAAAATGCACGATCCTACCGGATCTTGACGCGACAAAAGGAAGATCGGTTAGGCAAGTACGACACGGGGTATTCTACAAGTTGGGCAGAATGTGAGACATGGCTTATACCCCTATCACTATGTCGCGTTAGGTAGTTGATGAAATTATGGATACGGGGAAGAGCGATGTTGGGGTCGTGTTGTTGCCGTCCGGTGAACGGGCAATCATGCTCATGGGTACTGGTCACTGGTAGAATGTAATTGGCTCGTAGGGCCTCCGAACGCTGGCCTACCGTAATGGTGAATTTAGGCCAGAGAAACGGCAGATCTATGCTAGTTCTTTGCTATGCGTCCCCAATGCGCTGGCTGACGCGTTAAACACCCATCAGTGAAAGTCCCATAGTCCCGGATTCTCAGTATCTGACGCGTACTTGCCCATGTTGGACACCCTGCTGGAACGCCGTGCCTAAGCTGGGTTGTAAAACCATGGGGTGCCCCGGGTCGTTGGGTCCTGACCTAGTAGAAGTCTATGTCCGGCTTAAGGTACAGGACTACGGTTTATTGTAGTCCCAAGAAAGGTTGTGTTGGGGGCGGAAACTTATGGTGTGTGAGTCGGGGTTTCGTCGTAGGTATCGCAGTTTTCACAAAGGGAGCGTGCCCGACGCTTCTTCACGGATGACTCACTTCATAGGCTACTACCAACCGTCAATTCGTCTTGAAACCCAACCTCTGAAACCCATACAATTTAATGTACCGTACATGGATTTGTACAATTGCCGCGTGGATATCATTAGCAGAAGGCCACCTGTGGTCAGAGTGAGGCGGGTTGGCCTTTTCGGCTGTACCATGGACAGAATTAAGGAAAGCGTGGACTATCCTTACAGATACTCTGTTATAGGGTAAATTGTTCTCGTTTTTAGACTTGGGTATACATAGAGCAACCTAGAGCTGGAGTTATCACGGAGTGCTGATGAGAACTTTTTCTTTTCGACGAGGGACCAGTGAGGTCGGCACGCTAAGATTAGGGTTCAGTCAGCTTGTTGGTATCTACGCGCGCTGCATAGGCAGTGTCACAGTGGCCTTAACCTCACCTTGCTCCGGGGATCTCCAGAGGTATACTGCATGAGGCAAAGAGCCAAATAGACGGAATACTTACAATCTTAGACGCCAACAATACGCGCTACATCAATGCGTGGTCATATCCGTCCTTGGAGGGTCTGTCGAGTCGATCCTTACTCTCATCGTGGCTACGAACCCGTCGTAGATGGGTGGAGGATCGGCTGCGGAGAAGGCTACCAGTTCTCAGCTGACTGCGCAGATAATTCGTGGTATGTTGTCCAGCACAATGAGCCTGACAAGCATTTGGTACCGTCGCGTCAAAATGTGTAACTTGTAAAGATAATCACCAGCAGGGACTGATAGACTTAGCGTCAGGTCATCACGAAAAACCATTCCGCCGACGTATACGGTCGTTTAGCCCTACACAGGACAAATCCGAGGAGATAAATCATCTCCGGAACCCTGCGGGCAAACCCTCTGTGCTCTTAGTAAGGTTTTCCTGTTAGTAACATAGCGGAGGACGCCGTGACGATCTCAGGATTGCCGGCTACGCGGGTCCGAACTAGCTGACTATGCAAACTTTGAGTGGTGTCCGTGGTTCCAAAGAGGCGTGGGTTGAGTACGGTTCTGTAAAGTGCCTCTACTACGAACTGTCTAGATTAACCAGGAGGCCAACTACATTAAGCTATTATGACCATTGTATCAGGTTGGACCGCGCAGTCGTACAAGTTCTGCACTCCCGACTGGTTACAACCAAGATATGATACTAGCGCTGTATATGATCCCGAGAAAGGGTTATATCCTTACAGTAGGATCAGAAATCTTCCTCTGGGTCGAAATCCATTACTTA"
    #kmer_detector(Text,k)
    reverse_complement(sequence)

if __name__ == "__main__":
    main()
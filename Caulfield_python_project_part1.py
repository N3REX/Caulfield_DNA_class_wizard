#!/usr/bin/env python
# coding: utf-8

# Helpful code --------------------

standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}


# Sequence class --------------------


class seq:
    def __init__(self, name, organism, sequence, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    # define the info function
    def info(self):
        print(self.name)
        print(self.organism)
        print(self.sequence)
        print(self.type)

    # define the length function
    def length(self):
        lengthvar = len(self.sequence)
        print(lengthvar)

    # define the fasta_out function.
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + (self.name)
            + "_"
            + (self.organism)
            + "_"
            + (self.type)
            + "\n"
            + self.sequence
        )
        f.close()


# Protein class --------------------


class protein(seq):
    def __init__(self, name, organism, sequence, type, size):
        # add new attributes- i
        self.size = size  # in kDa

        # with super- no self, but need to add the arguments
        super().__init__(name, organism, sequence, type)

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "_Size="
            + self.size
            + "\n"
            + self.sequence
        )
        f.close()

    # methods of parent class are inherited
    # To the protein class, add a method called mol_weight, which returns the total molecular
    # weight of the protein sequence. The variable aa_mol_weights in the “helpful variables”
    # file should be helpful. This is a python dictionary of molecular weights for each amino
    # acid.
    def mol_weight(self):
        molsum = 0
        for amino in self.sequence:

            molsum += aa_mol_weights[amino]
        return molsum


# Nucleotide class --------------------


class nucleotide(seq):
    def __init__(self, name, organism, sequence, type):

        # with super- no self, but need to add the arguments
        super().__init__(name, organism, sequence, type)

    def gc_content(self):
        GCcount = 0
        for i in self.sequence:
            if i == "G" or i == "C":
                GCcount = GCcount + 1
        Bigcount = 0
        for i in self.sequence:
            Bigcount = Bigcount + 1
        Content = 100 * (GCcount / Bigcount)
        print(Content)

    # Laura's version was neater than mine
    # def gc_content(self):
    #    total = len(self.sequence)
    #    gc = self.sequence.count('G') + self.sequence.count('C')
    #    gc_percent = 100*gc/total
    #    print(gc_percent)


# DNA class --------------------


class DNA(nucleotide):
    def __init__(self, name, organism, sequence, type):

        # with super- no self, but need to add the arguments
        super().__init__(name, organism, sequence, type)

    def transcribe(self):
        Transcribed = self.sequence.replace("T", "U")  # Replaces Ts with Us
        return Transcribed

    # For the transcribe method
    # This site may be helpful https://www.geeksforgeeks.org/python-string-replace/

    # a. a method called six_frames that returns all 6 coding frames of the sequence (3
    # frames on the forward strand, 3 frames on the reverse complement strand)
    def six_frames(self):
        CodFram = []
        for i in range(3):
            CodFram.append(self.sequence[i:])
        nega_sequence = self.sequence[::-1]
        flip = ""
        for bp in nega_sequence:
            if bp == "A":
                flip += "T"
            elif bp == "T":
                flip += "A"
            elif bp == "G":
                flip += "C"
            elif bp == "C":
                flip += "G"
            else:
                flip += "N"
        for i in range(3):
            CodFram.append(flip[i:])
        return CodFram

    # b. a method called reverse_complement that returns the reverse complement of the
    # sequence
    def reverse_compliment(self):
        nega_sequence = self.sequence[::-1]
        # print(nega_sequence) #checking if the complement comes out right
        flip = ""
        for bp in nega_sequence:
            if bp == "A":
                flip += "T"
            elif bp == "T":
                flip += "A"
            elif bp == "G":
                flip += "C"
            elif bp == "C":
                flip += "G"
            else:
                flip += "N"
        return flip


# RNA class --------------------


class RNA(nucleotide):
    def __init__(self, name, organism, sequence, type):

        # with super- no self, but need to add the arguments
        super().__init__(name, organism, sequence, type)

    def start(self):
        print((self.sequence).find("AUG"))

    # For the start method
    # This site may be https://www.geeksforgeeks.org/python-string-find/

    def translate(self):
        # a. First, finds the start codon (AUG)
        start = (self.sequence).find("AUG")
        # b. Then, starting with AUG, breaks the sequence into 3 letter codons
        cones = [self.sequence[i : i + 3] for i in range(start, len(self.sequence), 3)]
        # c. Then translates those codons to amino acids
        brotein_bequence = ""
        for codon in cones:
            brotein_bequence += standard_code.get(codon, "-")
        # d. Returns the protein sequence of the amino acids
        return brotein_bequence


# Test --------------------

# Starting at the end of the notebook, assign the sequence:
# "CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA" to a variable of
# class DNA called uidA. This variable attributes should be name uidA, organism bacteria, and type DNA.
uidA = DNA(
    name="uidA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    organism="bacteria",
    type="DNA",
)

# 9. Use the fasta_out function to write the sequence and information for uidA DNA to a fasta file
uidA.fasta_out()

# 10. Use the six_frames function and reverse_complement functions to output the six coding
# frames and reverse complement of the uidA DNA sequence
uidA.reverse_compliment()
uidA.six_frames()

# 11. Using the transcribe function for the DNA class, transcribe the uidA DNA sequence to an RNA sequence.
uidA.transcribe()

# 12. Save this RNA sequence as a RNA class object called uidA_RNA with the same other
# attributes except the name should be uidA_RNA and the type should be RNA
uidARNASeq = uidA.transcribe()
uidA_RNA = RNA(
    name="uidA_RNA", sequence=str(uidARNASeq), organism="bacteria", type="RNA"
)

# 13. Use the fasta_out() function to write the RNA sequence and information for uidA to a fasta file
uidA_RNA.fasta_out()

# 14. Use the translate method on the RNA object uidA_RNA to translate the RNA sequence
uidA_RNA.translate()

# 15. Save this amino acid sequence as a protein class object called uidA_protein. Set the
# name as uidA_protein and the type as protein. You can set the size attribute as any
# value. Use the fasta_out() function to write this protein sequence and information to a new fasta file.
uidProt = uidA_RNA.translate()
uidA_protein = protein(
    name="uidA_protein",
    sequence=str(uidProt),
    organism="bacteria",
    type="protein",
    size="1337",
)
uidA_protein.fasta_out()

# 16. Use the method mol_weight to output the molecular weight of the amino acid sequence uidA_protein
uidA_protein.mol_weight()

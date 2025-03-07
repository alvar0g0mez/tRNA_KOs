"""
Define functions that might be needed in any of the other files of the project.
"""


def reverse_transcribe_nucleotide(nt: str):
    if nt == "A": return("T")
    elif nt == "C": return("G")
    elif nt == "G": return ("C")
    elif nt == "U": return("A")
    else:
        print("Provided letter is not a DNA nucleotide")


def transcribe_nucleotide(nt: str):
    if nt == "A": return("U")
    elif nt == "C": return("G")
    elif nt == "G": return ("C")
    elif nt == "T": return("A")
    else:
        print("Provided letter is not a DNA nucleotide")


def codon_to_anticodon(codon: str):
    anticodon = ""
    codon = codon[::-1]
    for nt in codon:
        new_nt = transcribe_nucleotide(nt)
        anticodon += new_nt
    return(anticodon)


def anticodon_to_codon(anticodon: str):
    codon = ""
    anticodon = anticodon[::-1]
    for nt in anticodon:
        new_nt = reverse_transcribe_nucleotide(nt)
        codon += new_nt
    return(codon)


































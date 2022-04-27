# accgcaatgtgggattcaacgacgcggaatagttaggtaaaa
# start: ATG
# stop: UAA, UAG, UGA

STOP = ['TAA', 'TAG', 'TGA']
START = 'ATG'
codonsMapping = {
    'UUU': 'Phe',
    'UUC': 'Phe',
    'UUA': 'Leu',
    'UUG': 'Leu',

    'UCU': 'Ser',
    'UCC': 'Ser',
    'UCA': 'Ser',
    'UCG': 'Ser',

    'UAU': 'Tyr',
    'UAC': 'Tyr',

    'UGU': 'Cis',
    'UGC': 'Cis',
    'UGG': 'Trp',

    # ====
    'CUU': 'Leu',
    'CUC': 'Leu',
    'CUA': 'Leu',
    'CUG': 'Leu',

    'CCU': 'Pro',
    'CCC': 'Pro',
    'CCA': 'Pro',
    'CCG': 'Pro',

    'CAU': 'His',
    'CAC': 'His',
    'CAA': 'Gln',
    'CAG': 'Gln',

    'CGU': 'Arg',
    'CGC': 'Arg',
    'CGA': 'Arg',
    'CGG': 'Arg',

    # ==
    'AUU': 'Ile',
    'AUC': 'Ile',
    'AUA': 'Ile',
    'AUG': 'Met',

    'ACU': 'Thr',
    'ACC': 'Thr',
    'ACA': 'Thr',
    'ACG': 'Thr',

    'AAU': 'Asn',
    'AAC': 'Asn',
    'AAA': 'Lys',
    'AAG': 'Lys',

    'AGU': 'Ser',
    'AGC': 'Ser',
    'AGA': 'Arg',
    'AGG': 'Arg',

    # ==
    'GUU': 'Val',
    'GUC': 'Val',
    'GUA': 'Val',
    'GUG': 'Val',

    'GCU': 'Ala',
    'GCC': 'Ala',
    'GCA': 'Ala',
    'GCG': 'Ala',

    'GAU': 'Asp',
    'GAC': 'Asp',
    'GAA': 'Glu',
    'GAG': 'Glu',

    'GGU': 'Gly',
    'GGC': 'Gly',
    'GGA': 'Gly',
    'GGG': 'Gly'
}


def getIndexStartCodon(sequence):
    result = sequence.find(START)
    while result != -1 and result % 3 != 0:
        result = sequence.find(START, result + 1)

    return result


def complement(codon):
    return codon.replace('A', 'U').replace('T', 'A').replace('C', '-').replace('G', 'C').replace('-', 'G')


def translateToProtein(label, sequence):
    sequence = sequence.upper()
    startCodonIndex = getIndexStartCodon(sequence)
    if startCodonIndex == -1:
        print("No genes found into the sequence!")
        return

    proteinGenesList = []
    proteinGene = ''
    codonIndex = startCodonIndex + 3
    while codonIndex < len(sequence):
        codon = sequence[codonIndex: codonIndex + 3]
        if codon in STOP:
            if proteinGene != '':
                proteinGenesList.append(proteinGene)
                proteinGene = ''

            sequence = sequence[codonIndex + 3:]
            startCodonIndex = getIndexStartCodon(sequence)
            if startCodonIndex == -1:
                break

            codonIndex = startCodonIndex + 3
            continue

        codon = complement(codon)
        proteinGene += ('-' if proteinGene != '' else '') + codonsMapping[codon]
        codonIndex += 3

    print("Printing found results for %s:" % label)
    if len(proteinGenesList) == 0:
        print("No genes found into the sequence!")
    else:
        n = 1
        for geneItem in proteinGenesList:
            print("%d# - %s" % (n, geneItem))
            n += 1
    print()


def main():
    translateToProtein('#1', 'accgcaatgtgggattcaacgacgcggaatagttaggtaaaa')
    translateToProtein('#2', 'accgcaatatgggattcaacgacgcggatgagttaggtaaaa')
    translateToProtein('#3', 'accgcaatatgggattcaacgacgcggatgagttatgttttt')
    translateToProtein('#4', 'atgccaatgaaatgagattcaacgacgcggaatatgaaggtataa')


if __name__ == "__main__":
    main()

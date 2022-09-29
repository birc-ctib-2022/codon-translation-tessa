"""Module for translating DNA to proteins via codons."""
#DICTIONARY reminder format for dictionary is d = { "key": "value", "key": "value"}
CODON_MAP = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
             'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
             'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
             'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
             'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
             'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
             'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
             'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
             'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
             'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
             'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
             'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}


def split_codons(dna: str) -> list[str] | None:
    #this function expects that the input will be a string called dna 
    #and that it will RETURN A LIST OF STRINGS 
    #REMINDER ()= string,[]= list 
    """Split a DNA string into a list of triplets.

    If the length of the string is a multiple of three, then this
    function splits the string into non-overlapping triplets.

    >>> split_codons('aaacccgggttt')
    ['aaa', 'ccc', 'ggg', 'ttt']

    If the string length is not a multiple of three, the function
    should return `None`. (There are better ways at reporting
    errors, but we will see those later).

    >>> split_codons("acgt") is None
    True

    """
    # FIXME: Implement the function
    codons= []
    #creates the empty list called codons that will be returned 
    if len(dna) % 3 != 0: 
        return None
        #if dna (the input string) is not divisible evenly by 3, return none 
    else: 
        for i in range(0,len(dna),3):
            #for i in range from 0 to length of string dna , taking 
            #steps of size 3 
            codons.append(dna[i:i+3])
            #list.append means take the empty list called list then the . calls 
            #the method append, and then inside the parenthesis are the arguments, 
            #in this case the arguments are take the string dna, 
            #add the first 3 indexes in that string as a new string 
            #then go back to the for loop that is the whole len(dna) and 
            #jump i 3 steps and repeat codons.append to create a new string () in the 
            #list called codons. 

    return codons


def translate_codons(codons: list[str]) -> list[str]:
    #input for function is a list [] of strings (), output should also be a LIST OF STRINGS
    """Translate a list of codons (triplets) into their corresponding
    amino acid sequence.

    >>> translate_codons(['TGT', 'TGC', 'TGA'])
    ['C', 'C', '*']

    The function must be able to handle both upper and lower case
    strings.

    >>> translate_codons(['tgt', 'tgc', 'tga'])
    ['C', 'C', '*']

    If the `codons` list contain anything that isn't a valid codon,
    i.e. not in the CODON_MAP when translated into upper case, the
    function should return `None`.

    >>> translate_codons(["acg", "ac", "gca"]) is None
    True

    """
    # FIXME: Implement the function
    aminoacids= []
    for i in codons:
        i=i.upper ()
        # call the method upper which is if input is lower case will convert it to upper case 
        if i not in CODON_MAP: 
            return None 
        else:
            aminoacids.append(CODON_MAP[i])
            #i is in brackets because we want to create a list 
    return aminoacids


def translate_dna(dna: str) -> str:
    #input is a string called dna, output is a string of amino acids, 
    #note that we need to have one string as our output but the in between functions both create lists of 
    #strings so i will need the .join method to joing the lists of strings into one string 
    """Translate a DNA string into its corresponding amino acid string.

    >>> translate_dna('TGTTGCTGA')
    'CC*'
    >>> translate_dna('tgttgctga')
    'CC*'

    If the sequence does not have a length that is a multiple of three, or if
    any of the triplets in it are not valid codons (when in uppercase), the function
    should return `None`.

    >>> translate_dna('tgtgctg') is None
    True

    """
    # FIXME: Implement the function
    if split_codons(dna) != None:
    #if the string dna is not length that is a multiple of 3 will return none 
        translate_codons(split_codons(dna))
        if translate_codons(split_codons(dna)) != None:
        #if any of the triplets is not a valid codon will return none 
            return "".join (translate_codons(split_codons(dna)))
            #the actual meat of the code .join method joins list of strings into one string 

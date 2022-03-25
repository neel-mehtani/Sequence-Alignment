# BWT.py
# HW1, Computational Genomics, Spring 2022
# andrewid: nmehtani

# WARNING: Do not change the file name, or the function signatures below.
# Autograder expects these names exactly.

def rle(s):
    """Run Length Encoder
    Args: s, string to be encoded
    Returns: RLE(s)
    """
    rle = ''

    i = 0
    #Iterate over every char in the string
    while i < len(s):
        char = s[i]
        j = i +1
        count = 1

        #Consider every char up to the last one
        if j < len(s) - 1:
          while s[j] == char:
            count += 1
            j += 1
          
          if count != 1:
            rle += char + char + f'{count}'
            i = j 
          else:
            rle += char
            i += 1

        #Last char
        else:
          rle += char
          i += 1
        

    return rle

def bwt_encode(s):
    """Burrows-Wheeler Transform
    Args: s, string, which must not contain '{' or '}'
    Returns: BWT(s), which contains '{' and '}'
    """
    s = "{" + s + "}"  # Add start and end of text marker
    
    table = sorted(s[i:] + s[:i] for i in range(len(s)))  # Table of rotations of string
    last_column = [row[-1:] for row in table]  # Last characters of each row
    
    return "".join(last_column)  # Convert list of characters into string

def bwt_decode(bwt):

    """Inverse Burrows-Wheeler Transform
    Args: bwt, BWT'ed string, which should contain '{' and '}'
    Returns: reconstructed original string s, must not contains '{' or '}'
    """

    table = [""] * len(bwt)  # Make empty table
    for i in range(len(bwt)):
        table = sorted(bwt[i] + table[i] for i in range(len(bwt)))  # Add a column of r
    s = [row for row in table if row.endswith("}")][0]  # Find the correct row (ending in ETX)
    
    return s.rstrip("}").strip("{")



def test_string(s):
    compressed = rle(s)
    bwt = bwt_encode(s)
    compressed_bwt = rle(bwt)
    reconstructed = bwt_decode(bwt)
    template = "{:25} ({:3d}) {}"
    print(template.format("original", len(s), s))
    print(template.format("bwt_enc(orig)", len(bwt), bwt))
    print(template.format("bwt_dec(bwt_enc(orig))", len(reconstructed), reconstructed))
    print(template.format("rle(orig)", len(compressed), compressed))
    print(template.format("rle(bwt_enc(orig))", len(compressed_bwt), compressed_bwt))
    print()
    print()

if __name__ == "__main__":
    # Add more of your own strings to explore for question (i)
    test_strings = ["WOOOOOHOOOOHOOOO!",
                    "scottytartanscottytartanscottytartanscottytartan"]
    for s in test_strings:
        test_string(s)

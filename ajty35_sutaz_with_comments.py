# Author: Miroslav Mocak
# Date: 21/FEB/2018
# Email: miroslav.mocak@ness.com

# Description: A single function accum:

# accum(“Python“)
# Výstup: “P-Yy-Ttt-Hhhh-Ooooo-Nnnnnn“

# accum(“IS“)
# Výstup: “I-Ss“

# accum(“cOmInG“)
# Výstup: “C-Oo-Mmm-Iiii-Nnnnn-Gggggg“
				
def accum(input):
    char = list(input) # convert input string to a list of letters
    out  = []          # define output list  
	
    for c in char:     # start of main logic: loop over every character in list
        d = c          # define temporary storage 
        for i in range(char.index(c)+1): 
            # append only the last string of the loop 
            # and capitalize first letter with the title method 
            if (i == char.index(c)): out.append(d.title()) 
            d += c # append character c to temporary storage

    print('-'.join(out)) # join all characters in the list with string separator '-'
		
accum("Python")	
accum("IS")	
accum("cOmInG")		
		
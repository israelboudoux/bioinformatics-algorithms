#############################
### DOTPLOT FUNCTIONS
#############################

# create matrix given dimensions: number of rows and columns filled with zeros
def create_mat(nrows, ncols):
    mat = []
    for i in range(nrows):
        mat.append([])
        for j in range(ncols):
            mat[i].append(0)
    return mat

# basic dotplot algorithm: fills with ones coincident characters
def dotplot(seq1, seq2):
    ''' Create a matrix based on the the input sequences and fill the cells that correspond to a match
    '''
    # ....
    return mat

# extended dotplot with window and stringency parameters
def extended_dotplot (seq1, seq2, window, stringency):
    mat = create_mat(len(seq1), len(seq2))
    start = int(window/2)
    for i in range(start,len(seq1)-start):
        for j in range(start, len(seq2)-start):
            matches = 0
            l = j - start
            for k in range(i-start, i+start+1):
                if seq1[k] == seq2[l]: matches += 1
                l += 1
                if matches >= stringency: mat[i][j] = 1
    return mat

# prints dotplot
def print_dotplot(mat, s1, s2):
    ''' Create a function to visualize the matrix
    Print each row in a different line
    if there is a match print the symbol "*" otherwise blankspace
    Print the symbols of seq2 as columns and seq1 as rows
    use the sys.stdout.write(...) to output
    '''
    import sys
    sys.stdout.write(" " + s2+"\n")
    # ....


def test_diagonal_length(mat, istart, jstart):
    # given the starting indices on the row and column
    # check along the diagonal that starts in istart and jstart
    # the longest sub-sequences of matches; return this value
    # ....


def test():
    s1 = "CGATATAGATT"
    s2 = "TATATAGTAT"
    mat1 = dotplot(s1, s2)
    print_dotplot(mat1, s1, s2)
    print(test_diagonal(mat1, 2, 3 ))
    #print("")
    #mat2 = extended_dotplot(s1, s2, 5, 4)
    #print_dotplot(mat2, s1, s2)


if __name__ == "__main__":
    test()

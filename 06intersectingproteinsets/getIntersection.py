#find the proteins that are both targeted for degradation and involved in division

#first clpXP-based, from Flynn et al. (2003).
filename = 'clpxpproteins.csv'

clpXPgenes = []

with open(filename, "r") as fo:
    for line in fo:
        clpXPgenes.append(line.split(';')[0].strip('?'))

#then clpAP-based. Convert to gene. Parsed from The N-degradome of Escherichia coli: limited proteolysis in vivo generates a large pool of proteins bearing N-degrons.
filename = 'clpAPproteins.txt'

clpAPgenes = []
tolower = lambda s: s[:1].lower() + s[1:] if s else ''

with open(filename, "r") as fo:
    for line in fo:
        clpAPgenes.append( tolower(line.split(' ')[0]).strip())

combo = set(clpAPgenes) | set(clpXPgenes)


#now intersect with known division genes. Obtained from Zhou, J. & Rudd, K. E. (2013).
filename = 'divisiongenes.txt'

division_genes = []

with open(filename, "r") as fo:
    for line in fo:
        division_genes.append(line.split('\t')[1].strip())

intersection = set(division_genes)  & set(combo)
print( intersection)
print(len(intersection))
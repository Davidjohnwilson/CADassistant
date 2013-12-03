clausestr = "[[2,3],[0,1,4]]"
polys = ["00000","111111","222222","333333","44444"]


clausenums = []   #clauses as numbers
clauselist = clausestr[1:-1].replace('[','|')
clauselist = clauselist.replace(']','|')
clauselist = clauselist.split('|') #split according to closing brackets
print(clauselist)
clauselist = [i for i in clauselist if i != '']
clauselist = [i for i in clauselist if i != ',']
print(clauselist)
for i in xrange(len(clauselist)):   #ignore final bracket
    indivclause = []
    indivclauselist = clauselist[i].split(',')
    print(indivclauselist)
    for j in xrange(len(indivclauselist)):
        print(indivclauselist[j])
        indivclause.append(int(indivclauselist[j]))
    clausenums.append(indivclause)
print(clausenums)

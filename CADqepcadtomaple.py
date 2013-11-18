import sys

arglist = sys.argv[1:]

#Ask user for the file to process
#inputfileanddir = raw_input('QEPCAD to Maple convertor \n ~~~~~~~~~~~~~~~~~~~~ \n Enter Input Filename (with directory if necessary): \n')

inputfileanddir = arglist[0]

#open file, extract lines and close file
f = open(inputfileanddir,'r')
flines = f.readlines()
f.close()

#compute number of lines in file
numlines = len(flines)

#initialise the output
outputstr =  'F,vars:='

	
#extract the name, variables, and polynomials
qname = flines[0]
qvars = flines[1]
qnumvars = flines[2]
qpolys = flines[3] 

mname = '#' + qname[1:-2] + '\n'

mvars = qvars[1:-2].split(',')
mvars = mvars[::-1]
mvars = ','.join(mvars)
mvars = '['+mvars+']'

mpolys = qpolys[:-2]
mpsplit = list(mpolys)
mpslen = len(mpsplit)
mpolynew = '['

k=0
while (k < mpslen):
    if mpsplit[k] == 'A' or mpsplit[k] == 'E':
        while mpsplit[k] != ')':
            k=k+1
        k=k+1
    elif mpsplit[k] == '[':
        mpolynew = mpolynew
        k=k+1
    elif mpsplit[k] == ']':
        mpolynew = mpolynew + ','
        k=k+3
    elif mpsplit[k] == '/' and mpsplit[k+1] == '\\':
        mpolynew = mpolynew
        k=k+2
    elif mpsplit[k] == '\\':
        mpolynew = mpolynew
        k=k+2
    elif mpsplit[k] == '=' and mpsplit[k+1] == ' ':
        if mpsplit[k+2] == '-':
            mpolynew = mpolynew + '+'
        else:
            mpolynew = mpolynew + '-'
        k = k+2
        #while (mpsplit[k] != ']'):
        #    mpolynew = mpolynew + mpsplit[k]
        #    k=k+1
    elif mpsplit[k] == '<' or mpsplit[k] == '>':
        if mpsplit == '=':
            k=k+1
        else:
            if mpsplit[k+2] == '(' and mpsplit[k+3] == '-':
                mpolynew = mpolynew + '+'
                k = k+1
            else:
                mpolynew = mpolynew + '-'
            k = k+2
            #while (mpsplit[k] != ']'):
            #    mpolynew = mpolynew + mpsplit[k]
            #    k=k+1
    elif mpsplit[k] == '~':
        k = k+1
    elif mpsplit[k] == '(' and (mpsplit[k+1] == 'E' or mpsplit[k+1] == 'A'):
        k = k+1


    
    elif mpsplit[k] == ' ':
        if len(mpolynew)>0 and (mpolynew[-1] == '-' or mpolynew[-1] == '+' or mpolynew[-1] == '*' or mpolynew[-1] == '/'): 
            k= k+1
        elif mpsplit[k+1] == '[':
            k = k+2
        elif mpsplit[k+1] == '+':
            mpolynew = mpolynew + '+'
            k=k+3
        elif mpsplit[k+1] == '-':
            mpolynew = mpolynew + '-'
            k=k+3
        elif mpsplit[k+1] == '=' or mpsplit[k+1] == '<' or mpsplit[k+1] == '>':
            k=k+1
        else:
            mpolynew = mpolynew + '*'
            k=k+1

    else:
        mpolynew = mpolynew + mpsplit[k]
        k=k+1


for k in xrange(len(mpolynew)-1):
    if mpolynew[k] == ',' and (mpolynew[k+1] == '*' or mpolynew[k+1] == ',' or mpolynew[k+1] == '/' or mpolynew[k+1] == '+' or mpolynew[k+1] == '-' ):
        mpolynew = mpolynew[:k+1] + ' ' + mpolynew[k+2:]

for k in xrange(len(mpolynew)-2):
    if mpolynew[k]=='=' and mpolynew[k+1]=='-':
        mpolynew = mpolynew[:k] + '  ' + mpolynew[k+2:]

for k in xrange(len(mpolynew)-1):
    if mpolynew[k]=='=':
        mpolynew = mpolynew[:k] + ' ' + mpolynew[k+1:]
        

klast = len(mpolynew)-1
while mpolynew[klast] == ' ' or mpolynew[klast] == ',':
    klast = klast - 1
if mpolynew[klast] == '*':
    klast = klast - 1
mpolynew = mpolynew[:klast+1]



mpolynew = mpolynew + ']'


minput = '\n' + mname + mpolynew + ',' + mvars

outputstr = outputstr + minput

outputstr = outputstr + ':\n'

featuresstring = """
NoOfPolysWithVar:=proc(F,v)
local tot,f:
tot := 0:
for f in F do
  if degree(f,v)>0 then
    tot := tot+1;
  end if:
end do:
return tot:
end proc:

NoOfMonosWithVar:=proc(F,v)
local tot,f,m:
tot := 0:
for f in F do
  for m in f do
    if degree(m,v)>0 then
      tot := tot+1;
    end if:
  end do:
end do:
return tot:
end proc:

NoOfPossMonosD:=proc(deg,vars)
local tot,n,l,T,v:
tot := 0:
n:=nops(vars):

l:=[seq(i,i=0..deg)]:
T:=combinat:-cartprod([seq(l,i=1..n)]):

while not T[finished] do
  v := T[nextvalue]():
  if (add(i,i in v)=deg) then
    tot := tot + 1:
  end if:
end do:

return tot:
end proc:

NoOfPossMonos:=proc(deg,vars)
local tot,i:
tot:=0:
for i from 0 to deg do
  tot := tot + NoOfPossMonosD(i,vars):
end do: 
return tot:
end proc:

FeatureComputation:=proc(F,vars)
  local Features:

  Features:=[]:

  #Number of Input Polys
  Features := [op(Features), nops(F)]:
  
  #Max Total Degree
  Features := [op(Features), max(seq(degree(f),f in F))]:

  #Max Degree in x0
  Features := [op(Features), max(seq(degree(f,x0),f in F))]:

  #Max Degree in x1
  Features := [op(Features), max(seq(degree(f,x1),f in F))]:

  #Max Degree in x2
  Features := [op(Features), max(seq(degree(f,x2),f in F))]:

  #Proportion of polynomials in which x0 occurs
  Features := [op(Features), NoOfPolysWithVar(F,x0)/nops(F)]:

  #Proportion of polynomials in which x1 occurs
  Features := [op(Features), NoOfPolysWithVar(F,x1)/nops(F)]:

  #Proportion of polynomials in which x2 occurs
  Features := [op(Features), NoOfPolysWithVar(F,x2)/nops(F)]:

  #Max Norm of the polynomials
  Features := [op(Features), max(seq(max({seq(abs(c),c in coeffs(expand(f),vars))}),f in F))]:

  #Proportion of monomials in which x0 occur
  Features := [op(Features), NoOfMonosWithVar(F,x0)/add(nops(f),f in F)]:

  #Proportion of monomials in which certain variables occur
  Features := [op(Features), NoOfMonosWithVar(F,x1)/add(nops(f),f in F)]:

  #Proportion of monomials in which certain variables occur
  Features := [op(Features), NoOfMonosWithVar(F,x2)/add(nops(f),f in F)]:

  return Features: 
  
end proc:


print(FeatureComputation(F,vars)):


"""

heuristicstring = """
####Heuristics

with(RegularChains):
read("ProjectionCAD.mm"):
with(ProjectionCAD):

interface(warnlevel=0):

print(VariableOrderingHeuristic( vars, F , heuristic=S, SeeAll=true ));
print(VariableOrderingHeuristic( vars, F , heuristic=N, SeeAll=true ));
print(VariableOrderingHeuristic( vars, F , heuristic=BrownBasic ));

"""

outfile = open(arglist[1], 'w')
outfile.write(outputstr)
outfile.write(featuresstring)
outfile.write(heuristicstring)
outfile.close()

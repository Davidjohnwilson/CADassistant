# CAD Assistant Program

# Have access to the command line arguments
import sys

# Store them in arglist
arglist = sys.argv[1:]

welcome = """
############################################################
#     CAD Assistant                                        #
#          by David Wilson (D.J.Wilson@bath.ac.uk)         #
#                                                          #
#     A tool for computing CADs with various technologies  #
#   and techniques.                                        #
#                                                          #
############################################################

"""

# Create CAD error exception
class CADError( Exception ): pass

# CadProblem Class
# Core class of a CAD problem consisting of it's name, polynomials
# and variables. Includes procedure to output a printable string 
# of the problem.
class CadProblem:
    name = ""
    polys = []
    variables = []

    def __init__(self,name,polys,variables):
        self.name = name
        self.polys = polys
        self.variables = variables

    def printCAD(self):
        printstr = ""
        printstr += "Name : \t" + self.name + "\n"
        printstr += "Polys: \t" + ", ".join(self.polys) + "\n"
        printstr += "Vars : \t" + ",".join(self.variables) + "\n"
        return printstr

# CadProblemMethod Class (extends CadProblem)
# A CAD problem along with choices of how to construct it, what
# it needs to be invariant with respect to, and whether it uses
# any subCAD techniques. Includes procedures to set these
# parameters as well as produce printable strings for the CAD
# acronym and the overall problem.
class CadProblemMethod(CadProblem):
    def __init__(self,name,polys,variables,constr,inva,subCAD):
        CadProblem.__init__(self,name,polys,variables)
        self.constr = constr
        self.inva   = inva
        self.subCAD = subCAD

    def setConstr(self,constr):
        self.constr = constr
    def setInva(self,inva):
        self.inva   = inva
    def setSubCAD(self,subCAD):
        self.subCAD = subCAD

    def printCADproblem(self):
        printstr = self.printCAD()
        printstr += "CAD : \t" + self.CADacronym()
        return printstr  + "\n"

    def CADacronym(self):
        #Simple function to return correct abbreviation
        CADstring = "CAD"
        if self.inva.lower() == "sign" or self.inva.lower() == "si":
            CADstring = "SI"+CADstring
        elif self.inva.lower() == "order" or self.inva.lower() == "oi":
            CADstring = "OI"+CADstring
        elif self.inva.lower() == "equational constraint" or self.inva.lower() == "ec":
            CADstring = "EC"+CADstring
        elif self.inva.lower() == "truth table" or self.inva.lower() == "tti":
            CADstring = "TTI"+CADstring

        if self.constr.lower() == "projection" or self.constr.lower() == "pl":
            CADstring = "PL-"+CADstring
        elif self.constr.lower() == "triangular sets" or self.constr.lower() == "regular chains" or self.constr.lower() == "t" or self.constr.lower() == "rc":
            CADstring = "RC-"+CADstring

        if self.subCAD.lower() == "manifold" or self.subCAD.lower() == "m":
            CADstring = "M-"+CADstring
        elif self.subCAD.lower() == "layered" or self.subCAD.lower() == "l":
            CADstring = "M-"+CADstring
        elif self.subCAD.lower() == "layered manifold" or self.subCAD.lower() == "lm":
            CADstring = "LM-"+CADstring

        return CADstring #subCAD+"-"+constr+"-"+inva+"CAD"

# CadTTICAD Class (extends CadProblemMethod)
# A CAD Problem with the above choices and lists describing the 
# clauses and equational constraints. Includes procedures to check
# validity of choices, set choices, and produce a printable string 
# describing the CAD problem.
class CadTTICAD(CadProblemMethod):
    def __init__(self,name,polys,variables,constr,inva,subCAD,clauses,eqcons):
        CadProblemMethod.__init__(self,name,polys,variables,constr,inva,subCAD)
        self.clauses = clauses
        self.eqcons = eqcons
        #self.checkClauses()

    #Check input is a valid set of clauses - shouldn't be used directly by user!
    def checkClausesInp(self,clauses):
         for i in clauses:
            for j in i:
                if j>len(self.polys):
                    raise CADError("Clause references non-existent poly")

    #Check the defined clauses are correct
    def checkClauses(self):
        self.checkClausesInp(self.clauses)

    def setClauses(self,clauses):
        self.checkClausesInp(clauses)
        self.clauses = clauses

    def checkEqConsInp(self,eqcons):
        if len(self.clauses) != len(self.eqcons):
            CADError("Number of clauses != Number of EqCons!")
        for i in xrange(len(self.clauses)):
            for j in eqcons[i]:
                if j not in self.clauses[i]:
                    CADError("Equational Constraint chosen is not part of clause!")

    def checkEqCons(self):
        self.checkEqConsInp(self.eqcons)

    def setEqCons(self,eqcons):
        self.eqcons = eqcons
        self.checkEqConsInp(eqcons)

    def printTTICADClauses(self):
        printstr = "Clauses:\t"
        for cl in self.clauses:
            printstr += "[" + ",".join(map(str,cl)) + "]\t"
        printstr += "\nEqCons: \t"
        for ec in self.eqcons:
            printstr += "[" + ",".join(map(str,ec)) + "]\t"
        return printstr + "\n"



def manualMethod():
    print(welcome)
    namestr = raw_input("Please enter a name for the CAD problem:")
    polystr = raw_input("Please enter the polynomials for your problem,\n within square brackets and separated by commas:")
    varstr  = raw_input("Please enter the variables for your problem,\n within square brackets and separated by commas:")
    # Separate the inputted polynomials
    polylist = polystr.split(',')
    polylist[0] = polylist[0][1:]     # Strip leading '['
    polylist[-1] = polylist[-1][:-1]  # Strip trailing '['
    # Separate the inputted variables
    varlist  = varstr.split(',')
    varlist[0] = varlist[0][1:]       # Strip leading '['
    varlist[-1] = varlist[-1][:-1]    # Strip trailing '['

    # Initialize CADProblemMethod object
    currCADproblem = CadProblemMethod(namestr,polylist,varlist,"","","")
    print("\nCurrent Problem:\n"+currCADproblem.printCADproblem())

    constrstr = raw_input("Please specify a construction method (default: projection&lifting):")
    if constrstr=="":
        constrstr = "PL"
    currCADproblem.setConstr(constrstr)

    invastr = raw_input("Please specify an invariance to build the CAD with respect to (default: sign-invariance):")
    if invastr=="":
        invastr = "SI"
    currCADproblem.setInva(invastr)

    subCADstr = raw_input("Please specify any subCAD techniques to use (default: none):")
    if subCADstr=="":
        subCADstr = ""
    currCADproblem.setSubCAD(subCADstr)

    print("\nCurrent Problem:\n"+currCADproblem.printCADproblem())

    print("Checking whether we are dealing with Equational Constraints or TTICAD...\n")
    if currCADproblem.inva.lower() == "eq" or currCADproblem.inva.lower() == "tti":
        print("Dealing with Equational constraints. Now identifying clauses.\n")
        # Cast object into a CadTTICAD object.
        currCADproblem.__class__ = CadTTICAD

        print("\nPolynomials are:")
        for j in xrange(len(currCADproblem.polys)):
                print("Poly " + str(j) + ": " + str(currCADproblem.polys[j]))


        clausestr = raw_input("Please enter in the clauses of the TTICAD (as a list of lists identifying the correct polynomials), or if dealing with equational constraints leave empty: ")
        if clausestr == "":
            currCADproblem.setClauses([list(xrange(len(currCADproblem.polys)))])
            print(currCADproblem.clauses)
        else:
            clausenums = []   #clauses as numbers
            clauselist = clausestr[1:-1].replace('[','|')
            clauselist = clauselist.replace(']','|')
            clauselist = clauselist.split('|') #split according to closing brackets
            clauselist = [i for i in clauselist if i != '']
            clauselist = [i for i in clauselist if i != ',']
            for i in xrange(len(clauselist)):   #ignore final bracket
                indivclause = []
                indivclauselist = clauselist[i].split(',')
                for j in xrange(len(indivclauselist)):
                    indivclause.append(int(indivclauselist[j]))
                clausenums.append(indivclause)
            currCADproblem.setClauses(clausenums)

        EClist = []
        for i in xrange(len(currCADproblem.clauses)):
            print("\nClause " + str(i) + ":")
            for j in currCADproblem.clauses[i]:
                print("Poly " + str(j) + ": " + str(currCADproblem.polys[j]))
            clauseEC = raw_input("Select the equational constraint for this clause: ")
            EClist.append([clauseEC])
        currCADproblem.setEqCons(EClist)

        print("\n" + currCADproblem.printTTICADClauses())


        #eqconsstr = raw_input("")
        #currCADproblem.setEqCons([[1]])
        #print(currCADproblem.printTTICADClauses())


    else:
        print("Dealing with non-equational constraint based CAD")


if __name__ == "__main__":
    if len(arglist)>0:
        if arglist[0]=='interactive':
            print "Interactive mode enabled\n"
    else:
        manualMethod()



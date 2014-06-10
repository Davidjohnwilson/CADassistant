# CAD Assistant Program

# Have access to the command line arguments
import sys

# Store them in arglist
arglist = sys.argv[1:]

welcome = """
############################################################
#     CADAssistant                                         #
#          by David Wilson (D.J.Wilson@bath.ac.uk)         #
#                                                          #
#     A tool for computing CADs with various technologies  #
#   and techniques.                                        #
#                                                          #
############################################################

"""

helpmessage = """

This is Version 1.0 of CADAssistant.

CADAssistant can be run in 'manual' or 'interactive' mode (default is the latter) by specifying either choice in the
command line arguments.

The name, polynomials, and variables of the problem should be input.

The polynomials should be given as a list within square brackets and separated by commas. Multiplication of both
constants and variables should be given using the asterisk symbol: '*'. Raising a variable to a power should be given
with the carat symbol: '^'.

Variables should be given as a comma-separated list within square brackets.

For example:
3D Example
[x^2+y^2+z^2-1, 4*x*y*z-1]
[x,y,z]

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

    def listOfPolys(self):
        liststr = "["
        for i in xrange(len(self.polys)):
            liststr += self.polys[i] + ", "
        liststr = liststr[:-2] + "]"
        return liststr

    def listOfVariables(self):
        liststr = "["
        for i in xrange(len(self.variables)):
            liststr += self.variables[i] + ", "
        liststr = liststr[:-2] + "]"
        return liststr

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
        printstr += "\nEqCons :\t"
        for ec in self.eqcons:
            printstr += "[" + ",".join(map(str,ec)) + "]\t\t"
        return printstr + "\n"

    def printTTICAD(self):
        printstr = self.printCADproblem()
        printstr += self.printTTICADClauses()
        return printstr + "\n"

#Initialising CAD
def initialisationMethod():
    print(welcome)
    namestr = raw_input("Please enter a name for the CAD problem (or 'help'):")

    if namestr.lower() == "help":
        print(helpmessage)
        namestr = raw_input("Please now enter the name for the CAD problem:")
    polystr = raw_input("Please enter the polynomials for your problem,\n within square brackets and separated by commas:")
    varstr  = raw_input("Please enter the variables for your problem,\n within square brackets and separated by commas:")
    # Separate the inputted polynomials
    polylist = polystr.split(',')
    polylist[0] = polylist[0][1:]     # Strip leading '['
    polylist[-1] = polylist[-1][:-1]  # Strip trailing '['
    for i in xrange(len(polylist)):
        while (polylist[i][0] == ' '):
            polylist[i] = polylist[i][1:] # strip leading ' '
    # Separate the inputted variables
    varlist  = varstr.split(',')
    varlist[0] = varlist[0][1:]       # Strip leading '['
    varlist[-1] = varlist[-1][:-1]    # Strip trailing '['

    # Initialize CADProblemMethod object
    currCADproblem = CadProblemMethod(namestr,polylist,varlist,"","","")
    print("\nCurrent Problem:\n"+currCADproblem.printCADproblem())

    return currCADproblem

def polyToQEPCAD(polynom):
    polynom = polynom.replace('*', ' ')
    polynom = polynom.replace('+', ' + ')
    polynom = polynom.replace('-', ' - ')
    return polynom


#manualMethod - let's the user choose how to construct their CAD
def manualMethod():
    currCADproblem = initialisationMethod()

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


        clausestr = raw_input("Please enter in the clauses of the TTICAD (as a list of lists identifying the correct\n"
                              " polynomials), or if dealing with equational constraints leave empty: ")
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

        print("\n\nSummary of TTICAD Problem:\n\n" + currCADproblem.printTTICAD())


        #eqconsstr = raw_input("")
        #currCADproblem.setEqCons([[1]])
        #print(currCADproblem.printTTICADClauses())


    else:
        print("Dealing with non-equational constraint based CAD")

#Interactive Method tries to discover the best formulation
def interactiveMethod():
    currCADproblem = initialisationMethod()

    # VARIABLE HEURISTIC CHOICE

    variableheuristic = raw_input("Do you wish to use a heuristic to choose a variable ordering? [Y/n] : ").lower()
    if variableheuristic == "n" or variableheuristic == "no":
        variableheuristic = False
    else:
        variableheuristic = True

    if variableheuristic:
        print("Please load Projection CAD to use any of the following heuristics.\n")

        manualheuristic = raw_input("Do you wish to pick the heuristic manually [m], or automatically [a]? [m/A] :").lower()
        if manualheuristic == "m" or manualheuristic == "manual" or manualheuristic == "manually":
            manualheuristic = True
        else:
            manualheuristic = False

        if manualheuristic:
            heuristicchoice = raw_input("Please select a heuristic to use from:\n Brown's Heuristic [B], sotd [S], ndrr [N], or 1-Layered Heuristic [LAY]\n (the default is Brown's Heuristic). [B/s/n/lay] :").lower()
            if heuristicchoice == "s" or heuristicchoice == "sotd":
                heuristicchoice = "S"
            elif heuristicchoice == "n" or heuristicchoice == "ndrr":
                heuristicchoice = "N"
            elif heuristicchoice == "lay" or heuristicchoice == "layered" or heuristicchoice == "l":
                heuristicchoice = "LayeredHeuristic"
            else:
                heuristicchoice = "BrownBasic"
        else:
            heuristicchoice = raw_input("Do you want a quick heuristic that may be less accurate [Q], or a slower heuristic that is more accurate [S]? [Q/s] :").lower()
            if heuristicchoice == "s" or heuristicchoice == "slow":
                heuristicchoice = "LayeredHeuristic"
            else:
                heuristicchoice = "BrownBasic"

        print("Please copy the following code into your Maple window:")
        heurstr = "S:=VariableOrderingHeuristic(["
        for i in xrange(len(currCADproblem.variables)):
            heurstr += currCADproblem.variables[i] + ", "
        heurstr = heurstr[:-2]+"], ["
        for i in xrange(len(currCADproblem.polys)):
            heurstr += currCADproblem.polys[i] + ", "
        heurstr = heurstr[:-2]+"]"
        heurstr += ",heuristic='"
        heurstr += heuristicchoice
        heurstr += "',SeeAll=false);\n"
        print(heurstr)

        heurvars = raw_input("Please paste the output from the above command: ")

         # Separate the inputted variables
        heurvars = heurvars.split(',')
        heurvars[0] = heurvars[0][1:]       # Strip leading '['
        heurvars[-1] = heurvars[-1][:-1]    # Strip trailing '['
        currCADproblem.variables = heurvars
        print("Variable ordering has been successfully changed.")


    # ONE LAYERED CAD CHOICE

    strictineq = raw_input("Does your problem involve only strict inequalities? [y/N] : ").lower()
    if strictineq == "y" or strictineq == "yes":
        strictineq = True
    else:
        strictineq = False

    if strictineq:
        print("You should use a 1-layered CAD. Please load Projection CAD and use the following input: ")
        strictstr = "C:=LCAD( ["
        for i in xrange(len(currCADproblem.polys)):
            strictstr += currCADproblem.polys[i] + ", "
        strictstr = strictstr[:-2] + "], 1, [" #remove final comma
        for i in xrange(len(currCADproblem.variables)):
            strictstr += currCADproblem.variables[i] + ", "
        strictstr = strictstr[:-2] + "]):\n"
        strictstr += "nops(C);"
        print(strictstr)

    # GROBNER PRECONDITIONING CHOICE

    grobnerpre = raw_input("Does your problem involve just a conjunction of equalities that you can use Grobner preconditioning for? [y/N] : ").lower()
    if grobnerpre=="y" or grobnerpre=="yes":
        grobnerpre = True
    else:
        grobnerpre = False

    if grobnerpre:
        tnoitest = raw_input("Do you wish to use TNoI to predict whether preconditioning will be useful? [Y/n] : ").lower()
        if tnoitest=="n" or tnoitest=="no":
            tnoitest = False
        else:
            tnoitest = True

        if tnoitest:
            print("Please run the following two lines of code: ")
            print("TNoI:=proc(F): add(nops(indets(f)),f in F): end proc:")
            print("TNoIdiff := TNoI(" + currCADproblem.listOfPolys() + ") - TNoI(Groebner[Basis](" + currCADproblem.listOfPolys() + ", plex(op(" + currCADproblem.listOfVariables() + "))));")
            grobner = raw_input("Is the value positive? [y/n] : ").lower()
            if grobner=="y" or grobner=="yes":
                grobner = True
            elif grobner=="n" or grobner=="no":
                grobner = False
        else:
            grobner = raw_input("Do you want to use Grobner preconditioning? [y/N] : ").lower()
            if grobner=="y" or grobner=="yes":
                grobner = True
            else:
                grobner = False

        if grobner:

            print("Please copy the following code into Maple to compute the Grobner basis:")

            print("Groebner[Basis](" + currCADproblem.listOfPolys() + ", plex(op(" + currCADproblem.listOfVariables() + ")));")

            polystr = raw_input("Please copy the output:")

            polylist = polystr.split(',')
            polylist[0] = polylist[0][1:]     # Strip leading '['
            polylist[-1] = polylist[-1][:-1]  # Strip trailing '['
            for i in xrange(len(polylist)):
                while (polylist[i][0] == ' '):
                    polylist[i] = polylist[i][1:] # strip leading ' '
            currCADproblem.polys = polylist

    # QUANTIFIER ELIMINATION CHOICE

    qeProblem = raw_input("Does your problem require quantifier elimination? [y/N] : ").lower()
    if qeProblem == "y":
        qeProblem = True
    else:
        qeProblem = False

    if qeProblem:
        print("You should use Qepcad or ... ")
        quantifiers = raw_input("Please input the quantifiers [in the format '(A x)(E y)']:")
        numquants = quantifiers.count(')')
        numfree = len(currCADproblem.variables)-numquants

        qestr = "[ " + currCADproblem.name + ' ]\n'
        qestr += "("
        for i in xrange(len(currCADproblem.variables)):
            qestr += currCADproblem.variables[i] + ", "
        qestr = qestr[:-2] + ")\n"
        qestr += str(numfree) + "\n"

        print("For each polynomial please enter a relational operator [=, <, >, <=, >=, =/=].")
        polyrelations = []
        for i in xrange(len(currCADproblem.polys)):
            relation = raw_input("Relation for " + currCADproblem.polys[i] + " : ")
            polyrelations.append(relation)

        print("For each consecutive pair of polynomials please enter a boolean relation [/\\, \\/, ==>, <==, <==>].")
        boolrelations = []
        for i in xrange(len(currCADproblem.polys)-1):
            relation = raw_input("Relation for [" + currCADproblem.polys[i] + "] ??? [" + currCADproblem.polys[i+1] + "] : ")
            boolrelations.append(relation)

        #TODO: Negation '~'

        qepcadstr = quantifiers + "[ "
        for i in xrange(len(currCADproblem.polys)-1):
            qepcadstr += "[ " + polyToQEPCAD(currCADproblem.polys[i]) + " " +  polyrelations[i] + " 0 ]"
            qepcadstr += " " + boolrelations[i] + " "
        if len(currCADproblem.polys)>1:
            qepcadstr += "[ " + polyToQEPCAD(currCADproblem.polys[-1]) + " " + polyrelations[-1] + " 0 ]"
        qepcadstr += " ]."

        qestr += qepcadstr

        print("Copy the following input into Qepcad:\n")
        print(qestr)





if __name__ == "__main__":
    if len(arglist)>0:
        if arglist[0].lower() == 'interactive':
            print "Interactive mode enabled.\n"
            interactiveMethod()
        else:
            print "Manual mode enabled.\n"
            manualMethod()
    else:
        print "No method specified. Default: interactive mode enabled.\n"
        interactiveMethod()

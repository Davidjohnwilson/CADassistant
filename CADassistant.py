#CAD Assistant Program

# Have access to the command line arguments
import sys

#Store them in arglist
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


#Create CAD error exception
class CADError( Exception ): pass

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

# newCAD = CadProblem("Parabola", ["a*x^2 + b*x + c"], ["x","a","b","c"])
# print(newCAD.printCAD())


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



class CadTTICAD(CadProblemMethod):
    def __init__(self,name,polys,variables,constr,inva,subCAD,clauses,eqcons):
        CadProblemMethod.__init__(self,name,polys,variables,constr,inva,subCAD)
        self.clauses = clauses
        self.eqcons = eqcons
        self.checkClauses()

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
        self.checkEqConsInp(eqcons)
        self.eqcons = eqcons

    def printTTICADClauses(self):
        printstr = "Clauses:\t"
        for cl in self.clauses:
            printstr += "[" + ",".join(map(str,cl)) + "]\t"
        printstr += "\nEqCons: \t"
        for ec in self.eqcons:
            printstr += "[" + ",".join(map(str,ec)) + "]\t"
        return printstr + "\n"


#TTICAD tests
Tticad = CadTTICAD("Parabola", ["a*x^2 + b*x + c", "sdasd","asdasd"], ["x","a","b","c"],"PL","SI","L",[[0,1],[2]],[[1],[]])
#print(Tticad.CADacronym())
#Tticad.checkClausesInp([[1,2]])
#Tticad.checkEqCons()
print(Tticad.printTTICADClauses())











def manualMethod():
    print(welcome)
    namestr = raw_input("Please enter a name for the CAD problem:")
    polystr = raw_input("Please enter the polynomials for your problem,\n within square brackets and separated by commas:")
    varstr  = raw_input("Please enter the variables for your problem,\n within square brackets and separated by commas:")
    polylist = polystr.split(',')
    polylist[0] = polylist[0][1:]
    polylist[-1] = polylist[-1][:-1]
    varlist  = varstr.split(',')
    varlist[0] = varlist[0][1:]
    varlist[-1] = varlist[-1][:-1]
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


if __name__ == "__main__":
    if len(arglist)>0:
        if arglist[0]=='interactive':
            print "Interactive mode enabled\n"
    else:
        print "Manual mode enabled\n"
        manualMethod()



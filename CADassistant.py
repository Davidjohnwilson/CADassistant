#CAD Assistant Program

# Have access to the command line arguments
import sys

#Store them in arglist
arglist = sys.argv[1:]

if __name__ == "__main__":
    if len(arglist)>0:
        if arglist[0]=='interactive':
            print "Interactive mode enabled\n"
    else:
        print "Manual mode enabled\n"



def CADname(constr,inva,subCAD):
    #Simple function to return correct abbreviation
    CADstring = "CAD"

    if inva.lower() == "sign" or inva.lower() == "si":
        CADstring = "SI"+CADstring
    elif inva.lower() == "order" or inva.lower() == "oi":
        CADstring = "OI"+CADstring
    elif inva.lower() == "equational constraint" or inva.lower() == "ec":
        CADstring = "EC"+CADstring
    elif inva.lower() == "truth table" or inva.lower() == "tti":
        CADstring = "TTI"+CADstring

    if constr.lower() == "projection" or constr.lower() == "pl":
        CADstring = "PL-"+CADstring
    elif constr.lower() == "triangular sets" or constr.lower() == "regular chains" or constr.lower() == "t" or constr.lower() == "rc":
        CADstring = "RC-"+CADstring

    if subCAD.lower() == "manifold" or subCAD.lower() == "m":
        CADstring = "M-"+CADstring
    elif subCAD.lower() == "layered" or subCAD.lower() == "l":
        CADstring = "M-"+CADstring
    elif subCAD.lower() == "layered manifold" or subCAD.lower() == "lm":
        CADstring = "LM-"+CADstring

    return CADstring #subCAD+"-"+constr+"-"+inva+"CAD"

print(CADname("PL","TTI","LM")+"\n")


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

newCAD = CadProblem("Parabola", ["a*x^2 + b*x + c"], ["x","a","b","c"])
print(newCAD.printCAD())


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

    def printCADMethod(self):
        printstr = self.printCAD()
        printstr += "CAD : \t" + CADname(self.constr,self.inva,self.subCAD)
        return printstr  + "\n"



newCADmethod = CadProblemMethod("Parabola", ["a*x^2 + b*x + c","a*x-1"], ["x","a","b","c"],"PL","TTI","LM")
print(newCADmethod.printCADMethod())

newCADmblank = CadProblemMethod("Parabola", ["a*x^2 + b*x + c","a*x-1"], ["x","a","b","c"],"","","")
print(newCADmblank.printCADMethod())
newCADmblank.setConstr("RC")
newCADmblank.setInva("TTI")
newCADmblank.setSubCAD("M")
print(newCADmblank.printCADMethod())

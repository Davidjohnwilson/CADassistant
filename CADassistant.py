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


def choosevariable(problem,variables):
    if len(variables) == 1:
        return variables


# Projection CAD Package (V3.9) 
# University of Bath
# Maintained by Matthew England (M.England@bath.ac.uk)
# 19th November 2013  
#print("This is V3.9 of the ProjectionCAD module from 19th November 2013, designed and tested for use in Maple 16 and 17."):

ProjectionCAD := module()

description "ProjectionCAD module: Calculating various cylindrical algebraic decompositions via projection and lifting.":

#NOTE:  Maybe chnage either the name or behaviour of global variable PolySets

# Telling Maple that the module is designed to be a Maple package
option package: 

export
	CADDist, 
	CADNormDist, 
	CADProjection,
	CADLifting,
	CADFull,
	CADGenerateStack,
	ECCAD,
	ECCADProjFactors,
	ECCADProjOp,
	ECCADFormulations,
	ECCADHeuristic,
	LCAD, 
	LCADRecursive, 
	LCADDisplay, 
	LTTICAD, 
	MCADLiftOverLowCAD,
	MCAD,
	MTTICAD, 
	LMCAD,
	LMTTICAD,
	NumCellsInPiecewiseCAD,
	NumCellsInCAD,
	TTICAD,
	TTICADDist,
	TTICADNormDist,
	TTICADResCADSet,
	TTICADResCAD,
	TTICADProjFactors,
	TTICADProjOp,
	TTICADQFFFormulations,
	TTICADFormulations,
	TTICADQFFHeuristic,
	TTICADHeuristic,
	VariableOrderings,
	VariableOrderingHeuristic,
	sotd,
	ndrr:

local
	PCAD_SortEqsRC,
	PCAD_LdCf,
	PCAD_MakeMonic,
	PCAD_MakeInteger,
	PCAD_RemoveConstantMultiples,
	PCAD_SFBasis,
	PCAD_ContSet,
	PCAD_PrimSet,
	PCAD_SetGCD,
	PCAD_SetFactors,
	PCAD_DoIntervalsIntersect,
	PCAD_RCBtoCube,
	PCAD_DoCubesIntersect,
	PCAD_IsIntervalInsideAnother,
	PCAD_IsCubeInsideAnother,
	PCAD_CCADmat,
	PCAD_CCADpsc,
	PCAD_CCADpscset,
	PCAD_CCADreductum,
	PCAD_CCADreductumset,
	PCAD_CCADProj1,
	PCAD_CCADProj2,
	PCAD_CCADProj,
	PCAD_CCADProjPolys,
	PCAD_McCallumProj,
	PCAD_McCcoeffset,
	PCAD_McCProjPolys,
	PCAD_IsCellZeroDim,
	PCAD_IsCellFullDim,
	PCAD_MinimalDelineatingPolynomial,
	PCAD_CanRCHaveSolutionAtSP,
	PCAD_CanRSHaveSolutionAtSP,
	PCAD_MakeCoprimeOverPoint,
	PCAD_WhichRCInSplitHasSolutionAtSP,
	PCAD_MakeSquareFreeOverPoint,
	PCAD_SPtoRootOf,
	PCAD_SortCadCell,
	PCAD_GenerateStack,
	PCAD_ProjCADLift,
	PCAD_LWRCADtoPWCAD,
	TTI_ECReducedProjectionOperator,
	TTI_TTIProjectionOperator,
	TTI_ECCAD,
	TTI_ECPP,
	TTI_TTIPP,
	TTI_TTICAD,
	TTI_UnivariateCase,
	TTI_LowerDimCAD,
	TTI_TTIGenerateStack,
	TTI_ResCAD,
	TTI_ResCADSet,
	TTI_ExclProj_phi,
	Formulation_VarOrdFree,
	Formulation_VarOrdRestricted,
	Formulations_RGS, 
	Formulations_Classification,
	Formulations_SetPartitions,
	Formulations_FunctionPartitions,
	Formulations_ECsForTTI,
	Formulations_ECCAD,
	Formulations_QFF,
	Formulations_TTICAD,
	Heuristics_sotdP,
	Heuristics_sotdL,
	Heuristics_ndrr,
	Heuristics_BHFeatures,
	Heuristics_PickSuggestions,
	Heuristics_RunOnPPs,
	Heuristics_FindPos,
	Heuristics_Display,
	Heuristics_ECCAD,
	Heuristics_TTICADQFF,
	Heuristics_TTICAD,
	Heuristics_VariableOrdering,
	Heuristics_GVOCADFull,
	Heuristics_GVOCADFullBlocks,
	Heuristics_Brown,
	Heuristics_BrownBlocks,
	LCAD_ProjCADLiftOneLayered,
	LCAD_ProjCADLiftNLayered,
	LCAD_IsCellNLayeredDim,
	LCAD_IsCellRepOneLayeredDim,
	LCAD_IsCellRepNLayeredDim,
	LCAD_GenerateSeparateProj,
	LCAD_LWRCADtoPWCAD,
	LCAD_Distribution,
	LCAD_NormDistribution,
	LCAD_Recursive,
	LCAD_RecursiveLR,
	LTTI_TTICADNLayered,
	LTTI_TTINLayeredGenerateStack,
	LTTI_NLayeredUnivariateCase,
	LTTI_NLayeredLowerDimCAD,
	LTTI_Dist,
	LTTI_NormDist,
	MCAD_LiftToManifold,
	MCAD_ConstructLowDimCAD,
	MCAD_ManifoldCAD,
	MCAD_MTTICAD,
	LMCAD_ConstructLowDimLCAD,
	LMCAD_LMCAD,
	LMTTICAD_LMTTICAD:

#####################################################################
###  CONTENTS
###  Section 1: User Commands = Just error checking
###  Section 2: Miscellaneous Commands
###  Section 3: Collins Projection
###  Section 4: McCallum projection
###  Section 5: Lifting
###  Section 6: ECCAD and TTICAD
###  Section 7: Formulations for CAD
###  Section 8: Heuristics
###  Section 9: SubCAD (David Wilson)
#####################################################################

#####################################################################
###  Section 1: User Commands 
#####################################################################

CADProjection := proc(in_F::{list, set}, vars::list, {method::symbol:='McCallum'}, $) :: 'set':
local F:
description "CADProjection: Computes a set of projection polynomials for use in constructing a CAD.",
			"Input: A list or set of polynomials and a list of variables. The polynomials are defined over the variables with coefficients in Q. Also, an optional input specifying the method to use.",
			"Output: A list of polynomials.":

F := convert(in_F, 'set'); 
if F::'set(list(polynom))' then
	error("invalid input: %1 expects the first argument to be a list or set containing only polynomials.", procname): 
fi:
if vars::'Not(list(symbol))' or vars=[] then
	error("invalid input: %1 expects the second argument to be a list of symbols.", procname): 
fi:
F := expand(F):
if map(X->coeffs(X, vars), F)::'Not(set(rational))' then
	error("invalid input: %1 expects the polynomials in the first argument to be defined over the ring of variables in the second argument, with coefficients in the rationals", procname):
fi: 
if method::'Not(identical(Collins, McCallum))' then
	error("invalid input: expected method = Collins, or McCallum, received: method = %1", method)
fi:
if not has([_passed],[_options][1]) then 
	WARNING("no method was specified, McCallum's projection operator will be used"):
fi: 
if method='McCallum' then 
	return(PCAD_McCProjPolys(F,vars)):
elif method='Collins' then 
	return(PCAD_CCADProjPolys(F,vars)):
fi:
end proc:

#####################################################################

CADLifting := proc(in_F::{list, set}, vars::list, {failure::symbol:='warn', finalCAD::symbol:='SI', method::symbol:='McCallum', output::symbol:=list, retcad::nonnegint:=0}, $) :: 'set':
local F,dim:
description "CADLifting: Computes a CAD in list format from a set of projection polynomials.",
			"Input: A set of projection polynomials and a list of variables. Optional arguments specifying the method and whether the final output should be order invariant.",
			"Output: A CAD in list format. Each entry is an index and a sample point (encoded by a regular chain and isolating cube).":

F := convert(in_F, 'set'); 
if F::'Not(set(polynom))' then
	error("invalid input: %1 expects the first argument to be a list or set containing only polynomials.", procname): 
fi:
if vars::'Not(list(symbol))' then
	error("invalid input: %1 expects the second argument to be either a symbol, or list of symbols.", procname): 
elif vars=[] then
	error("invalid input: %1 expects the second argument to be either a symbol, or list of symbols.", procname): 
fi:
if map(X->coeffs(X, vars), F)::'Not(set(rational))' then
	error("invalid input: %1 expects the polynomials in the first argument to be defined over the ring of variables in the second argument, with coefficients in the rationals", procname):
fi: 
F := expand(F):
if method::'Not(identical(Collins, McCallum))' then
	error("invalid input: expected method = Collins, or McCallum, received: method = %1", method)
fi:
if not has([_passed],[_options][3]) then 
	WARNING("no method was specified. Assuming that the projection used was McCallum"):
fi:
if finalCAD::'Not'('identical(SI, OI)') then
	error("invalid input: expected finalCAD = SI or OI, received: finalCAD = %1", finalCAD)
fi:
if method='Collins' and finalCAD='OI' then 
	error("the final CAD can only be guaranteed order invariant when using McCallum algorithm."):
fi:
dim:=nops(vars):
if retcad>dim then 
	error("invalid input: %1 expects the optional argument retcad to be a non-negative integer less than or equal to the number of variables specified in the second argument.", procname): 
fi:
if output::'Not(identical(list, listwithrep, piecewise, rootof))' then
	error("invalid input: expected output = list, or listwithrep, or piecewise, or rootof, received: output = %1", output)
fi:
if failure::'Not(identical(warn, err, giveFAIL))' then
	error("invalid input: expected failue = warn, or err, or giveFAIL, received: failure = %1", output):
fi:
if not has(packages(),RegularChains) then 
	WARNING("the CAD has sample points encoded using regular chains. We recommend you load the RegularChains package for working with the output"): 
fi:
PCAD_ProjCADLift( F, vars, method, finalCAD, retcad, output, failure ):
end proc:

#####################################################################

CADFull := proc(in_F::{list, set}, vars::list, {failure::symbol:='warn', finalCAD::symbol:='SI', method::symbol:='McCallum', output::symbol:=list, retcad::nonnegint:=0}, $) :: 'set':
local F, pset, dim, Failure:
description "CADFull: Computes a CAD in list format from a set of projection polynomials.",
			"Input: A list or set of polynomials and a list of variables. The polynomials are defined over the variables with coefficients in Q. Also, an optional input specifying the method to use and whether the final output should be order invariant.",
			"Output: A CAD in list format. Each entry is an index and a sample point (encoded by a regular chain and isolating cube).":

F := convert(in_F, 'set'); 
if F::'Not(set(polynom))' then
	error("invalid input: %1 expects the first argument to be a list or set containing only polynomials.", procname): 
fi:
if vars::'Not(list(symbol))' then
	error("invalid input: %1 expects the second argument to be either a symbol, or list of symbols.", procname): 
elif vars=[] then
	error("invalid input: %1 expects the second argument to be either a symbol, or list of symbols.", procname): 
fi:
F := expand(F):
if map(X->coeffs(X, vars), F)::'Not(set(rational))' then
	error("invalid input: %1 expects the polynomials in the first argument to be defined over the ring of variables in the second argument, with coefficients in the rationals", procname):
fi: 
if method::'Not(identical(Collins, McCallum))' then
	error("invalid input: expected method = Collins, or McCallum, received: method = %1", method):
fi:
if not has([_passed],[_options][3]) then 
	WARNING("no method was specified, McCallum's algorithm will be used"):
fi:
if finalCAD::'Not'('identical(SI, OI)') then
	error("invalid input: expected finalCAD = SI or OI, received: finalCAD = %1", finalCAD):
fi:
if method='Collins' and finalCAD='OI' then 
	error("ehe final CAD can only be guaranteed order invariant when using McCallum algorithm."):
fi:
dim:=nops(vars):
if retcad>dim then 
	error("invalid input: %1 expects the optional argument retcad to be a non-negative integer less than or equal to the number of variables specified in the second argument.", procname): 
fi:
if output::'Not(identical(list, listwithrep, piecewise, rootof))' then
	error("invalid input: expected output = list, or listwithrep, or piecewise, or rootof, received: output = %1", output):
fi:
if failure::'Not(identical(warn, err, giveFAIL))' then
	error("invalid input: expected failue = warn, or err, or giveFAIL, received: failure = %1", output):
fi:
if finalCAD='OI' and has([_passed], failure)=false then
	Failure := 'giveFAIL':
else
	Failure := failure:
fi:
if not has(packages(),RegularChains) then 
	WARNING("the CAD has sample points encoded using regular chains. We recommend you load the RegularChains package for working with the output"): 
fi:
if method='McCallum' then 
	pset := PCAD_McCProjPolys(F,vars):
	userinfo(2, {'ProjectionCAD', 'CADFull'}, "produced set of", nops(pset), "projection factors using the", method, "algorithm."):
	userinfo(4, {'ProjectionCAD', 'CADFull'}, "the projection factors are:", pset):
elif method='Collins' then 
	pset := PCAD_CCADProjPolys(F,vars):
	userinfo(2, {'ProjectionCAD', 'CADFull'}, "produced set of", nops(pset), "projection factors using the", method, "algorithm."):
	userinfo(4, {'ProjectionCAD', 'CADFull'}, "the projection factors are:", pset):
fi:
PCAD_ProjCADLift( pset, vars, method, finalCAD, retcad, output, Failure ):
end proc:

#####################################################################

NumCellsInPiecewiseCAD := proc( pwcad::piecewise, $ ) :: posint:
local i,size,num,nest:
description "NumCellsInPiecewiseCAD: Calculate the number of cells in a MapleCAD given in the piecewise output format.",
			"Input: A CAD in piecewise format. This could be output from either this module or the Regular Chains implementation.",
			"Output: The number of cells in the CAD.":

size := nops(pwcad)/2: 
num := 0:
for i from 1 to size do
	nest := op(2*i,pwcad): 
	if nest::piecewise then 
		num := num + NumCellsInPiecewiseCAD( nest ):
	else 
		num := num+1:
	fi:
od: 
num:
end proc:

#####################################################################

NumCellsInCAD := proc( cad )
description "PCAD_NuMCellsInCAD: Gives the number of cells in a CAD of given output format.":

if cad::piecewise then
	WARNING("this command is assuming the input is a CAD in piecewise output format"):
	return(NumCellsInPiecewiseCAD(cad)):
else
	WARNING("this command is assuming the input is a CAD in either list, listwithrep, rootof or cadcell output format"):
	return(nops(cad)):
fi:
end proc:

#####################################################################

CADGenerateStack := proc( cell::list, in_F::{list, set}, vars::list, {output::symbol:=list}, $ ) :: 'list':
local F, R, cellIndex, cellSP, SPrc, SPcube, retstack, alpha, cellN:
description "CADGenerateStack: Generate a stack over a cell given a set of polynomials. Uses TRDgenerate_stack but first error checks and preconditions polynomials so they are of the right format.",
			"Input: A cell from a CAD in list or listwithrep output format, a list or set of polynomials and a list of variables.  The polynomials should be defined over the variables and the cell be from a cad of dimension 1 less."
			"Output: A list of new cells forming a stack over the original cell.":

if vars::'Not(list(symbol))' then
	error("invalid input: %1 expects the second argument to be either a symbol, or list of symbols.", procname): 
elif vars=[] then
	error("invalid input: %1 expects the second argument to be either a symbol, or list of symbols.", procname): 
fi:
R := RegularChains:-PolynomialRing(vars):
F := convert(in_F, 'list'); 
if F::'Not(list(polynom))' then
	error("invalid input: %1 expects the first argument to be a list or set containing only polynomials.", procname): 
fi:
F := expand(F):
if map(X->coeffs(X, vars), F)::'Not(list(rational))' then
	error("invalid input: %1 expects the polynomials in the first argument to be defined over the ring of variables in the second argument, with coefficients in the rationals", procname):
fi: 
if output::'Not(identical(list, listwithrep))' then
	error("invalid input: expected output = list, or listwithrep, received: output = %1", output)
fi:
cellN := nops(cell):
if cellN=2 then 
	if output='listwithrep' then
		error("invalid input: first argument is not a cell from a cad in listwithrep format"):
	fi:
elif cellN=3 then
	if output='list' then
		error("invalid input: first argument is not a cell from a cad in list format. If it is 'listwithrep' then specify output."):
	fi:
else
	error("invalid input: first argument is not a cell from a cad in list or listwithrep format, (incorrect number of entries)"):
fi:
cellIndex := cell[1]:
cellSP := cell[-1]:
if not cellIndex::'list(integer)' or nops(cellSP)<>2 then 
	error("invalid input: first argument is not a cell from a cad in list or listwithrep format, (entries not correct type)"):
fi:
SPrc := cellSP[1]:
SPcube := cellSP[2]:
if SPrc['type']<>'regular_chain' or SPcube::'Not(list(list))' then 
	error("invalid input: first argument is not a cell from a cad in list or listwithrep format, (SP entries not correct type)"):
fi: 
alpha := PCAD_SPtoRootOf( SPrc,SPcube, RegularChains:-PolynomialRing( vars[2..nops(vars)] ) ): 
F := remove(X->evalb(is( subs(alpha,X)=0 )), F):
retstack := PCAD_GenerateStack(cell, F, R, 'output'):
retstack := map(X->op(1,X),retstack): 
end proc:

#####################################################################

ECCADProjOp := proc( E::set, A::set, vars::list(symbol), $ ) :: set:
local mvar, lvars, chk:
description "ECCADProjOp: Computes one application of the reduced Projection operator for equational constraints.",
			"Input: Sets E and A,  the ordered list of variables.",  
			"Output: The set of polynomials obtained by applying the reduced projection operator.":

mvar := vars[1]:
lvars := remove(X -> X=mvar, vars):
if not E::'set'( 'polynom'( 'rational', vars )) then
	error("invalid input: %1 expects the entries of the set in the first argument to be polynomials in the variables stated.", procname):
fi:
chk := map(X -> has(X, vars[1]), E):
if chk<>{true} then
	error("invalid input: %1 expects the polynomials in the first set to all contain the main variable.", procname):
fi:
if not A::'set'( 'polynom'( 'rational', vars )) then
	error("invalid input: %1 expects the entries of the set in the first argument to be polynomials in the variables stated.", procname):
fi:
TTI_ECReducedProjectionOperator(E, A, mvar, lvars):
end proc:

#####################################################################

TTICADProjOp := proc( E::table, A::table, t::posint, vars::list(symbol), $ ) :: set:
local i, mvar, lvars, indE, indA, Epols, Apols, chk:
description "TTICADProjOp: Computes one application of the TTICAD Projection operator.",
			"Input: Tables E and A, an integer t, the ordered list of variables.  The integr t is the number of QFFs in Phi.",  
			"The tables have entries for 1..t with each entry a set of polynomials.",
			"Output: The set of polynomials obtained by applying the TTI projection operator.":

mvar := vars[1]:
lvars := remove(X -> X=mvar, vars):
indE := [indices(E, 'nolist')]:
if nops(indE)<>t or convert(indE,'set')<>{seq(i, i=1..t)} then
	error("invalid input: %1 expects the first argument to be a table indexed from 1 to t=%2.", procname, t):
fi:
indA := [indices(A, 'nolist')]:
if nops(indA)<>t or convert(indA,'set')<>{seq(i, i=1..t)} then
	error("invalid input: %1 expects the second argument to be a table indexed from 1 to t=%2.", procname, t):
fi:
Epols := map(X -> op(X), {entries(E, 'nolist')}):
if not Epols::'set'( 'polynom'( 'rational', vars )) then
	error("invalid input: %1 expects the entries of the table in the first argument to be polynomials in the variables stated.", procname):
fi:
chk := map(X -> has(X, vars[1]), Epols):
if chk<>{true} then
	error("invalid input: %1 expects the polynomials in the first table to all contain the main variable.", procname):
fi:
Apols := map(X -> op(X), {entries(A, 'nolist')}):
if not Apols::'set'( 'polynom'( 'rational', vars )) then
	error("invalid input: %1 expects the entries of the table in the first argument to be polynomials in the variables stated.", procname):
fi:
TTI_TTIProjectionOperator(E, A, t, mvar, lvars):
end proc:

#####################################################################

ECCAD := proc( IN::list, vars::list(symbol), {failure::symbol:='giveFAIL', LiftAll::boolean:=false, output::symbol:=list, retcad::nonnegint:=0}, $ ) 
local f, gs, dim:
description "ECCAD: Compute the CAD using equational constraint at first level.",
			"Input: A list containing f, a polynomial representing an equational constraint and gs, a list of polynomials representing the other constraints.",
			"Also a list of ordered variables and optionally, an ouput choice."
			"Output: A CAD which is sign invariant for f and also for g on those cells where f=0.":

if nops(IN)<>2 then
	error("invalid input: %1 expects the first argument to be a list with first entry polynomial and second entry list of polynomials.", procname):
else
	f := IN[1]:
	gs := IN[2]:
fi:
dim := nops(vars):
if retcad>dim then 
	error("invalid input: %1 expects the optional argument retcad to be a non-negative integer less than or equal to the number of variables specified in the second argument.", procname): 
fi:
if not f::'polynom' ('rational', vars) then
	error("invalid input: %1 expects each list in the first argument to be a polynomial in the stated variables over the rationals.", procname):
fi:
if not gs::'list'( 'polynom' ('rational', vars) ) then
	error("invalid input: %1 expects the second argument to be a list of polynomials in the stated variables over the rationals.", procname):
fi:
if output::'Not(identical(list, listwithrep, piecewise, rootof))' then
	error("invalid input: expected output = list, or listwithrep, or piecewise, or rootof, received: output = %1", output)
fi:
if failure::'Not(identical(warn, err, giveFAIL))' then
	error("invalid input: expected failue = warn, or err, or giveFAIL, received: failure = %1", output):
fi:
TTI_ECCAD( [expand(f), expand(gs)], vars, failure, LiftAll, output, retcad):
end proc:

#####################################################################

ECCADProjFactors := proc( IN::list, vars::list(symbol), $ ) :: set:
local f, gs:
description "ECProjFactors: Compute the projection factors used for projection by ECCAD.",
			"Input: f, a polynomial representing an equational constraint and gs, a list of polynomials representing the other constraints.",
			"Also a list of ordered variables."
			"Output: The set of projection factors.":

if nops(IN)<>2 then
	error("invalid input: %1 expects the first argument to be a list with first entry polynomial and second entry list of polynomials.", procname):
else
	f := IN[1]:
	gs := IN[2]:
fi:
if not f::'polynom' ('rational', vars) then
	error("invalid input: %1 expects each list in the first argument to be a polynomial in the stated variables over the rationals.", procname):
fi:
if not gs::'list'( 'polynom' ('rational', vars) ) then
	error("invalid input: %1 expects the second argument to be a list of polynomials in the stated variables over the rationals.", procname):
fi:
TTI_ECPP( [expand(f), expand(gs)], vars):
end proc:

#####################################################################

TTICAD:=proc( PHI::list, vars::list(symbol), {failure::symbol:='giveFAIL', LiftAll::boolean:=false, output::symbol:=list, retcad::nonnegint:=0}, $ ) 
local i, t, phi, f, gs, dim, ECs, PHIX, phix: 
description "TTICAD: Compute the TTICAD for a list of QFFs.",
			"Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial or list of polynomials (equational constraint) and the second a list of polynomials (other constraints).",
			"Also a list of ordered variables and optionally, an ouput choice and retcad level."
			"Output: A truth-table invariant CAD for PHI or failure if the polynomials are not well-oriented.":

t:=nops(PHI):
# Error checking
for i from 1 to t do 
	phi := PHI[i]:
	if not phi::'list' then 
		error("invalid input: %1 expects the first argument to be a list of lists (each representing a QFF).", procname):
	fi:
	if nops(phi)<>2 then 
		error("invalid input: %1 expects each list in the first argument to consist of two components; but list %2 has %3 components", procname, i, nops(phi)):
	fi:
	gs := expand(phi[2]):
	if not gs::'list'( 'polynom' ('rational', vars) ) then
		error("invalid input: %1 expects each list in the first argument to consist of two parts, with the second a list of polynomials in the stated variables.  This is not the case for QFF %2.", procname, i):
	fi:
	ECs:=phi[1]:
	if ECs::'polynom' then
		f := expand(phi[1]):
		if not f::'polynom' ('rational', vars) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first (the equational constraints) being polynomial in the stated variables.  This is not the case for QFF %2.", procname, i):
		fi:
	elif ECs::'list' then
		if not ECs::'list'( 'polynom' ('rational', vars) ) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first a list of polynomials in the stated variables.  This is not the case for QFF %2.", procname, i):
		fi:
	else
		error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first representing the equational constraint(s) either a polynomial or a list of polynomials.  This is not the case for QFF %2.", procname, i):
	fi:
od:
dim := nops(vars):
if retcad>dim then 
	error("invalid input: %1 expects the optional argument retcad to be a non-negative integer less than or equal to the number of variables specified in the second argument.", procname): 
fi:
if output::'Not(identical(list, listwithrep, piecewise, rootof))' then
	error("invalid input: expected output = list, or listwithrep, or piecewise, or rootof, received: output = %1", output)
fi:
if failure::'Not(identical(warn, err, giveFAIL))' then
	error("invalid input: expected failue = warn, or err, or giveFAIL, received: failure = %1", output):
fi:
# Dealing with undesignated ECs
PHIX := []:
for i from 1 to t do 
	phi := PHI[i]:
	ECs := phi[1]:
	if ECs::'list' then
		if nops(ECs)=0 then
			phix := [false, phi[2]]:
		elif nops(ECs)=1 then
			phix := [op(phi[1]), phi[2]]:
		else
			userinfo(4, {'ProjectionCAD'},"the user did not designate an equational constraint in QFF %1.  Hence the QFF Heuristic command is being used:", i):
			#NOTE:  This is not efficient as we are calculating the projection set twice!
			phix := op(Heuristics_TTICADQFF( phi, vars, 'NS', ['vary',80], 0, false, true, [1,1] )):
		fi:
	else 
		phix := phi:
	fi:
	PHIX := [op(PHIX), phix]:
od:
TTI_TTICAD( expand(PHIX), vars, failure, LiftAll, output, retcad):
end proc:

#####################################################################

TTICADProjFactors:=proc( PHI::list, vars::list(symbol), $ ) :: set:
local i, t, phi, f, gs, PHIX, phix, ECs: 
description "TTICADProjFactors: Compute the projection factors used for TTICAD for a list of QFFs.",
			"Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial or list of polynomials (equational constraint) and the second a list of polynomials (from the other constraints).",
			"Also a list of ordered variables."
			"Output: The set of projection factors.":

t:=nops(PHI):
# Error Checking
for i from 1 to t do 
	phi := PHI[i]:
	if not phi::'list' then 
		error("invalid input: %1 expects the first argument to be a list of lists (each representing a QFF).", procname):
	fi:
	if nops(phi)<>2 then 
		error("invalid input: %1 expects each list in the first argument to consist of two components; but list %2 has %3 components", procname, i, nops(phi)):
	fi:
	ECs:=phi[1]:
	if ECs::'polynom' then
		f := expand(phi[1]):
		if not f::'polynom' ('rational', vars) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first (the equational constraints) being polynomial in the stated variables.  This is not the case for QFF %2.", procname, i):
		fi:
	elif ECs::'list' then
		if not ECs::'list'( 'polynom' ('rational', vars) ) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first a list of polynomials in the stated variables.  This is not the case for QFF %2.", procname, i):
		fi:
	else
		error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first representing the equational constraint(s) either a polynomial or a list of polynomials.  This is not the case for QFF %2.", procname, i):
	fi:
	gs := expand(phi[2]):
	if not gs::'list'( 'polynom' ('rational', vars) ) then
		error("invalid input: %1 expects each list in the first argument to consist of two parts, with the second a list of polynomials in the stated variables.  This is not the case for QFF %2.", procname, i):
	fi:
od:
# Dealing with undesignated ECs
PHIX := []:
for i from 1 to t do 
	phi := PHI[i]:
	ECs := phi[1]:
	if ECs::'list' then
		if nops(ECs)=0 then
			phix := [false, phi[2]]:
		elif nops(ECs)=1 then
			phix := [op(phi[1]), phi[2]]:
		else
			error("invalid input: %1 expects the equational constraints to have already been designated, but QFF %2 contains %3 constraints.", procname, i, nops(ECs)):
		fi:
	else 
		phix := phi:
	fi:
	PHIX := [op(PHIX), phix]:
od:
TTI_TTIPP( expand(PHIX), vars):
end proc:

#####################################################################

TTICADResCADSet:=proc( PHI::list, vars::list(symbol), $ ) :: list:
local i, t, phi, f, gs, PHIX, phix, ECs: 
description "TTICADResCADSet: Compute the ResCADSet.",
			"Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial or list of polynomials (equational constraint) and the second a list of polynomials (from the other constraints).",
			"Also a list of ordered variables."
			"Output: The ResCADSet.":

t:=nops(PHI):
for i from 1 to t do 
	phi := PHI[i]:
	if not phi::'list' then 
		error("invalid input: %1 expects the first argument to be a list of lists (each representing a QFF).", procname):
	fi:
	if nops(phi)<>2 then 
		error("invalid input: %1 expects each list in the first argument to consist of two parts; a polynomial (equational constraint) and a list of polynomials (other constraints)", procname):
	fi:
	ECs:=phi[1]:
	if ECs::'polynom' then
		f := expand(phi[1]):
		if not f::'polynom' ('rational', vars) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first (the equational constraints) being polynomial in the stated variables.  This is not the case for QFF %2.", procname, i):
		fi:
	elif ECs::'list' then
		if not ECs::'list'( 'polynom' ('rational', vars) ) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first a list of polynomials in the stated variables.  This is not the case for QFF %2.", procname, i):
		fi:
	else
		error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first representing the equational constraint(s) either a polynomial or a list of polynomials.  This is not the case for QFF %2.", procname, i):
	fi:
	gs := expand(phi[2]):
	if not gs::'list'( 'polynom' ('rational', vars) ) then
		error("invalid input: %1 expects each list in the first argument to consist of two parts, with the second a list of polynomials in the stated variables.  This is not the case for QFF %2.", procname, i):
	fi:
od:
# Dealing with undesignated ECs
PHIX := []:
for i from 1 to t do 
	phi := PHI[i]:
	ECs := phi[1]:
	if ECs::'list' then
		if nops(ECs)=0 then
			phix := [false, phi[2]]:
		elif nops(ECs)=1 then
			phix := [op(phi[1]), phi[2]]:
		else
			error("invalid input: %1 expects the equational constraints to have already been designated, but QFF %2 contains %3 constraints.", procname, i, nops(ECs)):
		fi:
	else 
		phix := phi:
	fi:
	PHIX := [op(PHIX), phix]:
od:
TTI_ResCADSet(PHIX, vars):
end proc:

#####################################################################

TTICADResCAD:=proc( PHI::list, vars::list(symbol), {output::symbol:=list}, $ ) :: list:
local mvar, lvars, i, t, phi, f, gs, ECs, PHIX, phix: 
description "TTICADResCAD: Compute a TTICAD for a list of QFFs using the RESCAD method.",
			"Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial or list of polynomials (equational constraint) and the second a list of polynomials (from the other constraints).",
			"Also a list of ordered variables and optionally, an ouput choice."
			"Output: Failure if the reduced projection set is not well oriented.  Otherwise a CAD for PHI which is guarenteed truth table invariant provided that no equational constraint is nullified.":

mvar:=vars[1]:
lvars:=remove(has,vars,mvar):
if nops(vars)>1 then
	WARNING("The output from TTICADResCAD is only guaranteed truth table invariant if no equational constraint is nullified by a point in %1", lvars):
fi:
t:=nops(PHI):
# Error Checking
for i from 1 to t do 
	phi := PHI[i]:
	if not phi::'list' then 
		error("invalid input: %1 expects the first argument to be a list of lists (each representing a QFF).", procname):
	fi:
	if nops(phi)<>2 then 
		error("invalid input: %1 expects each list in the first argument to consist of two parts; a polynomial (equational constraint) and a list of polynomials (other constraints)", procname):
	fi:
	ECs:=phi[1]:
	if ECs::'polynom' then
		f := expand(phi[1]):
		if not f::'polynom' ('rational', vars) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first (the equational constraints) being polynomial in the stated variables.  This is not the case for QFF %2.", procname, i):
		fi:
	elif ECs::'list' then
		if not ECs::'list'( 'polynom' ('rational', vars) ) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first a list of polynomials in the stated variables.  This is not the case for QFF %2.", procname, i):
		fi:
	else
		error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first representing the equational constraint(s) either a polynomial or a list of polynomials.  This is not the case for QFF %2.", procname, i):
	fi:
	gs := expand(phi[2]):
	if not gs::'list'( 'polynom' ('rational', vars) ) then
		error("invalid input: %1 expects each list in the first argument to consist of two parts, with the second a list of polynomials in the stated variables.  This is not the case for QFF %2.", procname, i):
	fi:
od:
if output::'Not(identical(list, listwithrep, piecewise, rootof))' then
	error("invalid input: expected output = list, or listwithrep, or piecewise, or rootof, received: output = %1", output)
fi:
# Dealing with undesignated ECs
PHIX := []:
for i from 1 to t do 
	phi := PHI[i]:
	ECs := phi[1]:
	if ECs::'list' then
		if nops(ECs)=0 then
			phix := [false, phi[2]]:
		elif nops(ECs)=1 then
			phix := [op(phi[1]), phi[2]]:
		else
			error("invalid input: %1 expects the equational constraints to have already been designated, but QFF %2 contains %3 constraints.", procname, i, nops(ECs)):
		fi:
	else 
		phix := phi:
	fi:
	PHIX := [op(PHIX), phix]:
od:
TTI_ResCAD(PHIX, vars, output):
end proc:

#####################################################################

ECCADFormulations := proc( IN::list, $) :: list:
local ecs, necs:
description "ECCADFormulations: Create the list of choices for entry into ECCAD.",
			"Input: Two lists of polynomials, the first assumed equational constraints.",
			"Ouput: A list of possible formulations for ECCAD.  Each list has a polynomial (designated EC) and a list of other constraints.":

if nops(IN)<>2 then
	error("invalid input: %1 expects the first argument to be a list with first entry polynomial and second entry list of polynomials."):
else
	ecs:=IN[1]:
	necs:=IN[2]:
fi:
if not ecs::'list(polynom(rational))' then
	error("invalid input: %1 expects each list in the first argument to contain polynomials over the rationals."):
fi:
if not necs::'list'( 'polynom' ('rational') ) then
	error("invalid input: %1 expects each list in the first argument to contain polynomials over the rationals."):
fi:
Formulations_ECCAD( [expand(ecs), expand(necs)] ):
end proc:

##################################################################### map(X-> [X[1],X[2..nops(X)] ], 

TTICADQFFFormulations := proc( IN::list(list(polynom)), {splitting::boolean:=true}, $ ) :: list:
local F, G:
description "TTICADQFFFormulations: Create the list of possible formulations for a TTICAD QFF.",
			"Input: Two lists of polynomials, the first assumed the equational constraints.  Optionally, the choice of whether to consider splitting.",
			"Ouput: The list of all possible formulations for TTICAD, by default including splitting.":

F := IN[1]:
G := IN[2]:
if splitting=false then
	return( map(X->[X],Formulations_ECCAD(F,G)) ):
else 
	return( Formulations_QFF( [expand(F), expand(G)] ) ):
fi:
end proc:

#####################################################################

TTICADFormulations := proc( LL::list(list(list(polynom))), {splitting::boolean:=true}, $ ) :: list:
description "TTICADFormulations: Create the list of possible formulations for TTICAD.",
			"Input: A list of lists each representing a QFF.  Each QFF has two lists of polynomials, the first assumed the equational constraints.  Optionally, the choice of whether to consider splitting.",
					"Note: it is assumed that the QFFs come from a disjunctive normal form.  I.e. trivial merging is not possible.",
			"Ouput: The list of all possible formulations for TTICAD, by default including splitting.":

Formulations_TTICAD(LL, splitting):
end proc:

#####################################################################

sotd := proc(IN, $) :: posint:
description "sotd: Calculate the sotd (sum of total degree) of a polynomial or list / set of polynomials.",
			"Input: Either a polynomial of list / set of polynomials.",
			"Ouput: The sotd of the polynomials or sum of the sotds for a list / set of polynomials.":

if IN::'set'('polynom') or IN::'list'('polynom') then 
	return( Heuristics_sotdL(IN) ):
elif IN::'polynom' then
	return( Heuristics_sotdP(IN) ):
else
	error("invalid input:  sotd expects its input to be a polynomial, or list/set of polynomials"):
fi:
end proc:

#####################################################################

ndrr := proc(INN, var::symbol, {method::symbol:='unstated', naivedeg::posint:=80, timeout:='unstated'}, $) :: posint:
local TimeO:
description "ndrr: Calculate the ndrr (number of distinct real roots) of a polynomial or list / set of polynomials with respect to a variable.  Note: if multivariate polynomials are included then they are ignored.",
			"Input: Either a polynomial of list / set of polynomials, and a variable.  Options for the method used and a timeout to apply.",
			"Ouput: The ndrr.":

if method::'Not(identical(full, naive, vary, unstated, SFnaive))' then
	error("invalid input: expected method = full, or naive, or vary, received: method = %1", method)
fi:
if not INN::'{polynom, list(polynom), set(polynom)}' then 
	error("invalid input:  ndrr expects its input to be a polynomial, or list/set of polynomials"):
fi:
if timeout='unstated' then
	TimeO := 0:
elif timeout::posint then
	TimeO := timeout:
else
	error("invalid input:  ndrr expects the optional argument timeout to be a positive integer."):
fi:
if method='unstated' then
	return( Heuristics_ndrr(INN, var, 'SFnaive', naivedeg, false, TimeO) ):
else 
	return( Heuristics_ndrr(INN, var, method, naivedeg, true, TimeO) ):
fi:
end proc:

#####################################################################

ECCADHeuristic := proc(IN::list, vars::list(symbol), {heuristic::symbol:='NS', ndrrOptions::list:=['vary',80], ndrrTO:='unstated', SeeAll::boolean:=false, Wratio::list(nonnegint):=[1,1]}, $) :: list(list):
local ecs, necs, ndrrOp, NTimeO:
description "ECHeuristic: Selects formulation(s) for ECCAD using heuristics.",
			"Input: Two lists of polynomials (the first assumed equational constraints) and a variable ordering.  Optionally a heuristic choice (S, N, NS or SN, W) for which the default is NS.  If W is picked then Wratio should also be provided.",
			"Ouput: A single ECCAD formulation or list of the ECCAD formulations suggested by the heuristicselected partitions.":

if ndrrTO='unstated' then
	NTimeO := 0:
elif ndrrTO::posint then
	NTimeO := ndrrTO:
else
	error("invalid input:  the optional argument ndrrTO should be a positive integer."):
fi:
ecs := IN[1]:
necs := IN[2]:
if nops(IN)<>2 then
	error("invalid input: %1 expects the first argument to be a list with first entry polynomial and second entry list of polynomials."):
else
	ecs:=IN[1]:
	necs:=IN[2]:
fi:
if not ecs::'list(polynom(rational))' then
	error("invalid input: %1 expects each list in the first argument to contain polynomials over the rationals."):
fi:
if not necs::'list'( 'polynom' ('rational') ) then
	error("invalid input: %1 expects each list in the first argument to contain polynomials over the rationals."):
fi:
if heuristic::'Not(identical(N, S, NS, SN, W))' then
	error("expected heuristic = S, or N, or NS, or SN, or W, received: heuristic = %1", heuristic):
fi:
if nops(Wratio)<>2 then 
	error("The optional argument Wratio should be a list containing two nonnegative integers."):
fi:
if nops(ndrrOptions)>2 then
	error("The optional argument ndrrOptions should be a list with one or two entries"):
elif ndrrOptions[1]::'Not(identical(full, naive, vary))' then
	error("The optional argument ndrrOptions should be a list with first entry an acceptable ndrr method"):
elif nops(ndrr)=1 then
	ndrrOp:=[ndrrOptions[1], 80]:
elif not ndrrOptions[2]::posint then
		error("The optional argument ndrrOptions should be a list.  If the list contains a second entry then this should be a positive integer indicating the degree that the ndrr method changes."):
else
	ndrrOp:=ndrrOptions:
fi:
Heuristics_ECCAD( ecs, necs, vars, heuristic, ndrrOp, NTimeO, SeeAll, Wratio ):
end proc:

#####################################################################

TTICADQFFHeuristic := proc( IN::list(list(polynom)), vars::list(symbol), {heuristic::symbol:='NS', ndrrOptions::list:=['vary',80], ndrrTO:='unstated', SeeAll::boolean:=false, splitting::boolean:=true, Wratio::list(nonnegint):=[1,1]}, $) :: list(list):
local F, G, ndrrOp, NTimeO:
description "TTICADQFFHeuristic: Selects formulation(s) for a TTICADQFF using heuristics.",
			"Input: Two lists of polynomials (the first assumed equational constraints) and a variable ordering.  Optionally a heuristic choice (S, N, NS or SN) for which the default is NS and the choice to consider splittings (for which the default is true).",
			"Ouput: A single TTICAD formulation or list of the TTICAD formulations suggested by the heuristicselected partitions.":

if ndrrTO='unstated' then
	NTimeO := 0:
elif ndrrTO::posint then
	NTimeO := ndrrTO:
else
	error("invalid input:  the optional argument ndrrTO should be a positive integer."):
fi:
F := IN[1]:
G := IN[2]:
if heuristic::'Not(identical(N, S, NS, SN, W))' then
	error("expected heuristic = S, or N, or NS, or SN, or W, received: heuristic = %1", heuristic):
fi:
if nops(Wratio)<>2 then 
	error("The optional argument Wratio should be a list containing two nonnegative integers."):
fi:
if nops(ndrrOptions)>2 then
	error("The optional argument ndrrOptions should be a list with one or two entries"):
elif ndrrOptions[1]::'Not(identical(full, naive, vary))' then
	error("The optional argument ndrrOptions should be a list with first entry an acceptable ndrr method"):
elif nops(ndrr)=1 then
	ndrrOp:=[ndrrOptions[1], 80]:
elif not ndrrOptions[2]::posint then
		error("The optional argument ndrrOptions should be a list.  If the list contains a second entry then this should be a positive integer indicating the degree that the ndrr method changes."):
else
	ndrrOp:=ndrrOptions:
fi:
Heuristics_TTICADQFF( [expand(F), expand(G)], vars, heuristic, ndrrOp, NTimeO, SeeAll, splitting, Wratio ):
end proc:

#####################################################################

TTICADHeuristic := proc(LL::list(list(list(polynom))), vars::list(symbol), {modular::boolean:=false, heuristic::symbol:='NS', ndrrOptions::list:=['vary',80], ndrrTO:='unstated', SeeAll::boolean:=false, splitting::boolean:=true, Wratio::list(nonnegint):=[1,1]}, $) :: list(list):
local ndrrOp, NTimeO:
description "TTICADHeuristic: Selects formulation(s) for a sequence of QFFs for TTICAD using heuristics.",
			"Input: A list of lists each representing a QFF.  Each QFF has two lists of polynomials, the first assumed the equational constraints.",
					"Optionally a heuristic choice (S, N, NS or SN) for which the default is NS and the choice to consider splittings (for which the default is true).",
					"Note: it is assumed that the QFFs come from a disjunctive normal form.  I.e. trivial merging is not possible.",
			"Ouput: A single TTICAD formulation or list of the TTICAD formulations suggested by the heuristicselected partitions.":

if ndrrTO='unstated' then
	NTimeO := 0:
elif ndrrTO::posint then
	NTimeO := ndrrTO:
else
	error("invalid input:  the optional argument ndrrTO should be a positive integer."):
fi:
if heuristic::'Not(identical(N, S, NS, SN, W))' then
	error("expected heuristic = S, or N, or NS, or SN, or W, received: heuristic = %1", heuristic):
fi:
if nops(Wratio)<>2 then 
	error("The optional argument Wratio should be a list containing two nonnegative integers."):
fi:
if modular=true and SeeAll=true then 
	WARNING("SeeAll=true is specified but this is repressed when using the modular method.  Suggest running TTICADQFFHeuristic on the individual caluses."):
fi:
if nops(ndrrOptions)>2 then
	error("The optional argument ndrrOptions should be a list with one or two entries"):
elif ndrrOptions[1]::'Not(identical(full, naive, vary))' then
	error("The optional argument ndrrOptions should be a list with first entry an acceptable ndrr method"):
elif nops(ndrr)=1 then
	ndrrOp:=[ndrrOptions[1], 80]:
elif not ndrrOptions[2]::posint then
		error("The optional argument ndrrOptions should be a list.  If the list contains a second entry then this should be a positive integer indicating the degree that the ndrr method changes."):
else
	ndrrOp:=ndrrOptions:
fi:
Heuristics_TTICAD( LL, vars, modular, heuristic, ndrrOp, NTimeO, SeeAll, splitting, Wratio ):
end proc:

#####################################################################

VariableOrderings := proc( vars::list, $ ) :: list:
description "VariableOrderings: Finds all the possible acceptable variable orderings.",
			"Input: Either a list of variables (if there are no restrictions) or a list of lists of variables (if the variables must be in set blocks, e.g. for quantifier elimination.",
			"Ouput: A list of all posible acceptable variable orderings (each a list of variables).":

if vars::'list'(symbol) then
	return( Formulation_VarOrdFree(vars) ) :
elif vars::'list'('list'('symbol')) then
	return( Formulation_VarOrdRestricted(vars) ) :
else
	error("The input should either be a list of variables or a list of lists of variables.  Variables should be of type symbol."):
fi:
end proc:

#####################################################################

VariableOrderingHeuristic := proc( VarsList::list, IN, {algorithm::symbol:='CADFull', greedy::boolean:=false, heuristic::symbol:='NS', method::symbol:='McCallum', ndrrOptions::list:=['vary',80], ndrrTO:='unstated', SeeAll::boolean:=false, Wratio::list(nonnegint):=[1,1] }, $ ) :: list:
local vars, i, INN, f, gs, t, phi, alg, ndrrOp, NTimeO:
description "VariableOrderingHeuristic: Pick a variable ordering based on a heuristic.",
			"Input: The input for a CAD algorithm, a variable list (possibly in blocks) and an algorithm selection.",
			"Ouput: Suggested variable ordering(s).":

if ndrrTO='unstated' then
	NTimeO := 0:
elif ndrrTO::posint then
	NTimeO := ndrrTO:
else
	error("invalid input:  the optional argument ndrrTO should be a positive integer."):
fi:
if VarsList::'list(symbol)' then
	vars:=VarsList:
elif VarsList::'list(list(symbol))' then
	vars:=map(X->op(X),VarsList):
else 
	error("invalid input: %1 expects the first argument to be either a list of symbols, or list of lists of symbols.", procname): 
fi:
if heuristic::'Not(identical(N, S, NS, SN, W, BrownFull, BrownBasic))' then
	error("expected heuristic = S, or N, or NS, or SN, or W, or BrownFull, or BrownBasic, received: heuristic = %1", heuristic):
fi:
if (heuristic='BrownFull' or heuristic='BrownBasic') and algorithm<>'CADFull' then 
	error("The Brown heuristic is currently only implemented for algorithm=CADFull."):
fi:
if (heuristic='BrownFull' or heuristic='BrownBasic') and SeeAll=true then 
	error("The SeeAll option is repressed when using a greedy method such as Brown's heuristic."):
fi:
if algorithm=CADFull then
	if IN::'Not({list,set})' then
		error("invalid input: algorithm is set to %1 but the first argument is not admissible as a first argument for CADProjection", algorithm):
	fi:
	INN := convert(IN, 'set');
	if INN::'set(list(polynom))' then
		error("invalid input: %1 expects the first argument to be a list or set containing only polynomials.", algorithm): 
	fi:
	INN := expand(INN):
	if map(X->coeffs(X, vars), INN)::'Not(set(rational))' then
		error("invalid input: %1 expects the polynomials in the first argument to be defined over the ring of variables in the second argument, with coefficients in the rationals", algorithm):
	fi: 
	if method::'Not(identical(Collins, McCallum))' then
		error("invalid input: expected method = Collins, or McCallum, received: method = %1", method)
	fi:
	alg:='algorithm_CADFull':
elif algorithm=ECCAD then
	if IN::'Not(list)' then
		error("invalid input: algorithm is set to %1 but the first argument is not admissible as a first argument for ECCADProjFactors", algorithm):
	fi:
	if nops(IN)<>2 then
		error("invalid input: %1 expects the first argument to be a list with first entry polynomial and second entry list of polynomials.", algorithm):
	else
		f:=IN[1]:
		gs:=IN[2]:
	fi:
	if not f::'polynom' ('rational', vars) then
		error("invalid input: %1 expects each list in the first argument to be a polynomial in the stated variables over the rationals.", algorithm):
	fi:
	if not gs::'list'( 'polynom' ('rational', vars) ) then
		error("invalid input: %1 expects the second argument to be a list of polynomials in the stated variables over the rationals.", algorithm):
	fi:
	alg:='algorithm_ECCAD':
	INN := IN:
elif algorithm=TTICAD then
	if IN::'Not(list)' then
		error("invalid input: algorithm is set to %1 but the first argument is not admissible as a first argument for TTICADProjFactors", algorithm):
	fi:
	t:=nops(IN):
	for i from 1 to t do 
		phi := IN[i]:
		if not phi::'list' then 
			error("invalid input: %1 expects the first argument to be a list of lists (each representing a QFF).", algorithm):
		fi:
		if nops(phi)<>2 then 
			error("invalid input: %1 expects each list in the first argument to consist of two components; but list %2 has %3 components", algorithm, i, nops(phi)):
		fi:
		f := expand(phi[1]):
		gs := expand(phi[2]):
		if not f::'polynom' ('rational', vars) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first a polynomial in the stated variables.  This is not the case for QFF %2.", algorithm, i):
		fi:
		if not gs::'list'( 'polynom' ('rational', vars) ) then
			error("invalid input: %1 expects each list in the first argument to consist of two parts, with the second a list of polynomials in the stated variables.  This is not the case for QFF %2.", algorithm, i):
		fi:
	od:
	alg:='algorithm_TTICAD':
	INN := IN:
else
	error("expected algorithm = CADFull, or ECCAD, or TTICAD, received: algorithm = %1", algorithm):
fi:
if nops(Wratio)<>2 then 
	error("The optional argument Wratio should be a list containing two positive integers."):
fi:
if greedy=true and SeeAll=true then 
	WARNING("SeeAll=true is specified but this is repressed when using the greedy method.  Suggest running TTICADQFFHeuristic on the individual caluses."):
fi:
if nops(ndrrOptions)>2 then
	error("The optional argument ndrrOptions should be a list with one or two entries"):
elif ndrrOptions[1]::'Not(identical(full, naive, vary))' then
	error("The optional argument ndrrOptions should be a list with first entry an acceptable ndrr method"):
elif nops(ndrr)=1 then
	ndrrOp:=[ndrrOptions[1], 80]:
elif not ndrrOptions[2]::posint then
		error("The optional argument ndrrOptions should be a list.  If the list contains a second entry then this should be a positive integer indicating the degree that the ndrr method changes."):
else
	ndrrOp:=ndrrOptions:
fi:
Heuristics_VariableOrdering( VarsList, INN, alg, greedy, heuristic, method, ndrrOp, NTimeO, SeeAll, Wratio ):
end proc:

#####################################################################

LCAD:=proc(in_F::{list, set}, n::posint, vars::list, {failure::symbol:='warn', finalCAD::symbol:='SI', method::symbol:='McCallum', output::symbol:=list, retcad::nonnegint:=0}, $) :: 'set':
    local F, pset, dim, tmpCAD:
    description "LCAD: Computes an n-layered CAD in list format from a set of projection polynomials.",
                "Input: A list or set of polynomials, and a list of variables. The polynomials are defined over the variables with coefficients in Q. Also, an optional input specifying the method to use and whether the final output should be order invariant.",
                "Output: An n-layered CAD in list format. Each entry is an index and a sample point (encoded by a regular chain and isolating cube).":

    F := convert(in_F, 'set');
    if F::'Not(set(polynom))' then
        error("invalid input: %1 expects the first argument to be a list or set containing only polynomials.", procname):
    fi:
    if vars::'Not(list(symbol))' then
        error("invalid input: %1 expects the second argument to be either a symbol, or list of symbols.", procname):
    elif vars=[] then
        error("invalid input: %1 expects the second argument to be either a symbol, or list of symbols.", procname):
    fi:
    F := expand(F):
    if map(X->coeffs(X, vars), F)::'Not(set(rational))' then
        error("invalid input: %1 expects the polynomials in the first argument to be defined over the ring of variables in the second argument, with coefficients in the rationals", procname):
    fi:
    if method::'Not(identical(Collins, McCallum))' then
        error("invalid input: expected method = Collins, or McCallum, received: method = %1", method):
    fi:
    if not has([_passed],[_options][3]) then
        WARNING("no method was specified, McCallum's algorithm will be used"):
    fi:
    if finalCAD::'Not'('identical(SI, OI)') then
        error("invalid input: expected finalCAD = SI or OI, received: finalCAD = %1", finalCAD):
    fi:
    if method='Collins' and finalCAD='OI' then
        error("ehe final CAD can only be guaranteed order invariant when using McCallum algorithm."):
    fi:
    dim:=nops(vars):
    if retcad>dim then
        error("invalid input: %1 expects the optional argument retcad to be a non-negative integer less than or equal to the number of variables specified in the second argument.", procname):
    fi:
    if output::'Not(identical(list, listwithrep, piecewise, rootof))' then
        error("invalid input: expected output = list, or listwithrep, or piecewise, or rootof, received: output = %1", output):
    fi:
    if failure::'Not(identical(warn, err, giveFAIL))' then
        error("invalid input: expected failue = warn, or err, or giveFAIL, received: failure = %1", output):
    fi:
    if not has(packages(),RegularChains) then
        WARNING("the CAD has sample points encoded using regular chains. We recommend you load the RegularChains package for working with the output"):
    fi:
#    if not has(packages(),'ProjectionCAD') then
#        WARNING("the LayeredCAD package requires the ProjectionCAD module (V2.0+) which has not been loaded. Please load the package and rerun the command"):
#    fi:
#    if (kernelopts('opaquemodules')=true) then
#        WARNING("the subCAD package requires kernelopts(opaquemodules) to be false (to allow access to procedures within the RegularChains module). This will now be changed."):
#        kernelopts('opaquemodules'=false):
#    fi:
    if method='McCallum' then
        pset := PCAD_McCProjPolys(F,vars):
        userinfo(3, 'CADFull', "produced set of", nops(pset), "projection polynomials using the", method, "algorithm."):
    elif method='Collins' then
        pset := PCAD_CCADProjPolys(F,vars):
        userinfo(3, 'CADFull', "produced set of", nops(pset), "projection polynomials using the", method, "algorithm."):
    fi:
    if output='piecewise' then
        if n=1 then
            tmpCAD:=LCAD_ProjCADLiftOneLayered(pset,vars,method,finalCAD,retcad, 'listwithrep',failure):
        else
            tmpCAD:=LCAD_ProjCADLiftNLayered(pset,n,vars,method,finalCAD,retcad, 'listwithrep',failure):
        fi:
        return  LCAD_LWRCADtoPWCAD(tmpCAD):
    fi:
 
    if n=1 then
        return LCAD_ProjCADLiftOneLayered(pset,vars,method,finalCAD,retcad, output,failure):
    else
        return LCAD_ProjCADLiftNLayered(pset,n,vars,method,finalCAD,retcad, output,failure):
    fi:
end proc:

#####################################################################

CADDist := proc(F :: {list, set}, vars :: list, {method::symbol:='McCallum'}, $) :: 'list':
    description "CADDist: Generates a list with the number of cells of each dimension in a full CAD using McCallum projection.",
                "Input: A list or set of polynomials and a list of variables.",
                "Output: A list of the number of cells of each dimension.":
    if method::'Not(identical(Collins, McCallum))' then
        error("invalid input: expected method = Collins, or McCallum, received: method = %1", method):
    fi:
    return LCAD_Distribution(F,vars,method):
end proc:

#####################################################################

CADNormDist :=  proc(F :: {list, set}, vars :: list, {method::symbol:='McCallum'}, $) :: 'list':
    description "CADNormDist: Generates a list with the normalised number of cells of each dimension in a full CAD using McCallum projection.",
                "Input: A list or set of polynomials and a list of variables.",
                "Output: A list of the normalised number of cells of each dimension.":

    if method::'Not(identical(Collins, McCallum))' then
        error("invalid input: expected method = Collins, or McCallum, received: method = %1", method):
    fi:
    return LCAD_NormDistribution(F,vars,method):
end proc:

#####################################################################

LCADRecursive := proc(F::{list,set},vars::list,C::list,LD::list,{method::symbol:='McCallum',output::symbol:='list', resetproj::boolean:=true},$)
    description "LCADRecursive: Computes a recursive layered CAD in list or listwithrep format. ",
                "Input: A set of polynomials, list of ordered variables, list of cells a previous recursive call terminated at (possibly empty), a layered CAD in list format from a previous recursive call (possibly empty). Also, an optional parameter stating whether to reset the projection polynomials.",
                "Output: A layered CAD in list format of one more layer than inputted (or 1 layer if an empty CAD is given), an unevaluated recursive call to produce the following layer (using the 'value' call).":
    if output = 'list' then
        return LCAD_Recursive(F,vars,C,LD,method,output,resetproj):
    else
        return LCAD_RecursiveLR(F,vars,C,LD,method,output,resetproj):
    end if:
    # NOTE: NEEDS ERROR CHECK ON OUTPUT
end proc:

#####################################################################

LCADDisplay:=proc(LCAD::list) :: 'piecewise':
    description "CADDisplay: Converts a layered CAD in listwithrep format into a piecewise output (with placeholders for truncated branches).",
                "Input: A layered CAD in listwithrep format.",
                "Output: A piecewise construct for the layered CAD.":

    if LCAD::'Not(list)' then 
        error("invalid input: %1 expects a list as input", procname):
    fi:
    if nops(LCAD)=0 then
        error("invalid input: %1 expects a list as input", procname):
    fi:
    if nops(LCAD[1]) <> 3 then
        error("invalid input: %1 expects a CAD in listwithrep format", procname):
    fi:
    return LCAD_LWRCADtoPWCAD(LCAD):
end proc:

#####################################################################

LTTICAD:=proc( PHI::list, N, vars::list(symbol), {output::symbol:=list}, $ ) :: list:
    local i, t, phi, f, gs:
    description "LTTICAD: Compute the N-layered TTICAD for a list of QFFs.",
                "Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial (equational constraint) and the second a list of polynomials (from the other constraints).",
                "Also a list of ordered variables and optionally, an ouput choice."
                "Output: A truth-table invariant N-layered CAD for PHI or failure if the polynomials are not well-oriented.":

    t:=nops(PHI):
    for i from 1 to t do
        phi := PHI[i]:
        if not phi::'list' then
            error("invalid input: %1 expects the first argument to be a list of lists (each representing a clause).", procname):
        fi:
        if nops(phi)<>2 then
            error("invalid input: %1 expects each list in the first argument to consist of two components; but list %2 has %3 components", procname, i, nops(phi)):
        fi:
        f := expand(phi[1]):
        gs := expand(phi[2]):
        if not f::'polynom' ('rational', vars) then
            error("invalid input: %1 expects each list in the first argument to consist of two parts, with the first a polynomial in the stated variables.  This is not the case for clause %2.", procname, i):
        fi:
        if not gs::'list'( 'polynom' ('rational', vars) ) then
            error("invalid input: %1 expects each list in the first argument to consist of two parts, with the second a list of polynomials in the stated variables.  This is not the case for clause %2.", procname, i):
        fi:
    od:
    if output::'Not(identical(list, listwithrep, piecewise, rootof))' then
        error("invalid input: expected output = list, or listwithrep, or piecewise, or rootof, received: output = %1", output)
    fi:
    LTTI_TTICADNLayered( expand(PHI),N, vars, output):
end proc:

#####################################################################

TTICADDist := proc(PHI :: {list, set}, vars :: list) :: 'list':
    description "TTICADDist: Generates a list with the number of cells of each dimension in a full TTICAD using McCallum projection.",
                "Input: A list or set of polynomials and a list of variables.",
                "Output: A list of the number of cells of each dimension.":
    return LTTI_Dist(PHI,vars):
end proc:
 
#####################################################################

TTICADNormDist := proc(PHI :: {list, set}, vars :: list) :: 'list':
    description "TTICADNormDist: Generates a list with the normalised number of cells of each dimension in a full TTICAD using McCallum projection.",
                "Input: A list or set of polynomials and a list of variables.",
                "Output: A list of the normalised number of cells of each dimension.":

    return LTTI_NormDist(PHI,vars):
end proc:

#####################################################################

MCADLiftOverLowCAD:=proc(lowCAD::list,equCon,vars::list,{output::symbol:=list},$)::list:
    description "MCADLiftOverLowCAD: Lifts over an (n-1)-dimensional CAD to produce a Manifold CAD.",
                    "This is currently restricted to the case where the equational constraint being lifted with respect to has positive degree in the variable being lifted.",
                "Input: An (n-1)-dimensional CAD, an equational constraint to lift with respect to, and a list of variables.",
                "Output: A manifold sub-CAD which is a decomposition of the manifold defined by the equational constraint provided.":

    if degree(equCon,vars[1])<1 then
        ERROR("currently a Manifold sub-CAD can only be constructed for equational constraints with positive degree in %1.", vars[1]);
    fi:
    if nops(lowCAD)>0 and nops(lowCAD[1][1]) <> (nops(vars) -1) then
        ERROR("the lower-dimensional CAD provided is of dimension %1 and not of dimension %2 as required by the variables provided.", nops(lowCAD[1][1]), (nops(vars)-1)):
    fi:
    return MCAD_LiftToManifold(lowCAD,equCon,vars,output):
end proc:

#####################################################################

MCAD:=proc(T::list,vars::list)::list:
    description "MCAD: Produces a Manifold sub-CAD. Currently restricted to the case where the equational constraint being lifted with respect to has positive degree in the variable being lifted.",
                "Input: A list consisiting of an equational constraint to lift with respect to, and a list of non-equational constraints, and a list of variables.",
                "Output: A manifold sub-CAD which is a decomposition of the manifold defined by the equational constraint provided.":

    if degree(T[1],vars[1])<1 then
        ERROR("currently a Manifold CAD can only be constructed for equational constraints with positive degree in %1.", vars[1]);
    fi:
    return MCAD_ManifoldCAD(T,vars):
end proc:

#####################################################################

MTTICAD:=proc(PHI::list,vars::list)::list:
    local i:
    description "MTTICAD: Produces a Manifold TTICAD. Currently restricted to the case where the equational constraints being lifted with respect to have positive degree in the variable being lifted.",
                "Input: A list consisiting of QFFs, each with an equational constraint to lift with respect to, and a list of non-equational constraints; and a list of variables.",
                "Output: A manifold sub-TTICAD which is a decomposition of the manifold defined by the equational constraints provided.":

    for i from 1 to nops(PHI) do
        if degree(PHI[i][1],vars[1])<1 then
            ERROR("currently a Manifold CAD can only be constructed for equational constraints with positive degree in %1.", vars[1]);
        fi:
    od:
    return MCAD_MTTICAD(PHI,vars):
end proc:

#####################################################################

LMCAD:=proc(T::list,n::posint,vars::list)::list:
    description "LMCAD: Produces a n-layered manifold sub-CAD. Currently restricted to the case where the equational constraints being lifted with respect to have positive degree in the variable being lifted.",
                "Input: A list consisiting of an equational constraint to lift with respect to, and a list of non-equational constraints; a positive integer between 1 and nops(vars)+1 corresponding to the number of layers required; and a list of variables.",
                "Output: An n-layered manifold sub-CAD which is a decomposition of the manifold defined by the equational constraints provided into the top n layers.":

    if degree(T[1],vars[1])<1 then
        ERROR("currently a Layered Manifold CAD can only be constructed for equational constraints with positive degree in %1.", vars[1]);
    fi:

    if n < 1 or (n-1) > nops(vars) then
        ERROR("cannot make a %1 layered manifold CAD with %2 variables", n, nops(vars)):
    fi:
return LMCAD_LMCAD(T,n,vars):
end proc:

#####################################################################

LMTTICAD:=proc(PHI::list, n::posint, vars::list)::list:
    local i:
    description "LMTTICAD: Produces a Layered Manifold TTICAD. Currently restricted to the case where the equational constraints being lifted with respect to have positive degree in the variable being lifted.",
                "Input: A list consisiting of clauses, each with an equational constraint to lift with respect to, and a list of non-equational constraints; a positive integer between 1 and nops(vars)+1; and a list of variables.",
                "Output: An n-layered manifold sub-TTICAD which is a decomposition of the manifold defined by the equational constraints provided producing the top n layers.":

    for i from 1 to nops(PHI) do
        if degree(PHI[i][1],vars[1])<1 then
            ERROR("currently a Manifold CAD can only be constructed for equational constraints with positive degree in %1.", vars[1]);
        fi:
    end do:
    if n < 1 or (n-1) > nops(vars) then
        ERROR("cannot make a %1 layered manifold CAD with %2 variables", n, nops(vars)):
    fi:
    return LMTTICAD_LMTTICAD(PHI,n,vars):
end proc:

#####################################################################
#####################################################################
#####################################################################
### LOCAL COMMANDS
#####################################################################
#####################################################################
#####################################################################

#####################################################################
### Section 2: Miscellaneous useful things
#####################################################################

PCAD_SortEqsRC := proc( eqs, R, $ ) :: 'list':
local vars,dim,i,E:
description "PCAD_SortEqsRC: Sort a list of equations into ascending main variable. Error if main variables are not distinct.":

vars:=ListTools:-Reverse(R['variables']):
dim:=nops(vars):
E:=table():
for i from 1 to nops(vars) do
	E[i]:=select(X->RegularChains:-MainVariable(X,R)=vars[i], eqs): 
	if nops(E[i])>1 then 
		error("the equations inputted do not have distinct main variable. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
	fi:
od:
[seq( op(E[i]), i=1..dim)]
end proc: 

#####################################################################

PCAD_LdCf := proc( poly::polynom, var::symbol, $ ) :: 'polynom':
description "PCAD_LdCf: Computes the leading coefficient of a polynomial with respect to a variable":

coeff( poly, var, degree(poly,var)):
end proc:

#####################################################################

PCAD_MakeMonic:=proc( poly::polynom, RorVar, $ ) :: polynom:
local i, vars, R, pol: 
description "PCAD_MakeMonic: Divide polynomial by its leading coefficient (in the base field).":

if RorVar::'table' and RorVar['type']='polynomial_ring' then
	vars:=RorVar['variables']:
	R:=RorVar:
elif RorVar::'list'('symbol') then
	vars:=RorVar:
	R:=RegularChains:-PolynomialRing(vars):
else
	error("the second argument should either be a polynomial ring or a list of variables. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
if indets(poly,'name') minus convert(vars,set) <> {} then
	error("the polynomial has variables other than those specified. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
if poly::'constant' then 
	return(poly)
else 
	pol:=poly:
	i:=1:
fi:
while i<nops(vars)+1 do
	pol:=RegularChains:-Initial(pol,R):
	if pol::'constant' then 
		return( expand(poly/pol) ):
	else
		i:=i+1: 
	fi:
od:
error("leading coefficient could not be identified. This is an unexpected error message - please report to M.England@bath.ac.uk"):
end proc:

#####################################################################

PCAD_MakeInteger:=proc( poly::polynom, $ ) :: polynom:
local bottom:
description "PCAD_MakeInteger: If the denominator is a constant then multiply up.":

bottom:=denom(poly):
if bottom::'constant' then 
	return( expand(bottom*poly) ): 
else 
	return(poly) 
fi:
end proc:

#####################################################################

PCAD_RemoveConstantMultiples:=proc( pols, $ ) :: set;
local polys, vars, RET:
description "PCAD_RemoveConstantMultiples: Removes polynomials from a set that are constant multiples of polynomials previously appearing in the set.":

if pols::list then 
	polys:=convert(pols,'set'): 
elif pols::set then 
	polys:=pols:
else 
	error("first input should be a list or set of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
if polys={} then 
	return({}): 
fi:
vars:=convert(indets(polys,'name'),list):
RET := map(X->PCAD_MakeMonic(X, vars), polys):
RET := map(X->PCAD_MakeInteger(X), RET): 
end proc:

#####################################################################

PCAD_SFBasis:=proc( pols, var::symbol, $ ) :: set;
local polys, i, j, decomp, head, body, RET:
description "PCAD_SFBasis: Construct a square free basis of a set of polynomials, with respect to a variable":

if pols::list then 
	polys:=convert(pols,set): 
elif pols::set 
	then polys:=pols:
else 
	error("first input should be a list or set of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
RET:={}:
for i from 1 to nops(polys) do 
	decomp:=sqrfree(op(i,polys),var): 
	head:=op(1,decomp): 
	body:=op(2,decomp):
	if type(head,constant)=false then 
		RET:={op(RET),op(1,decomp)}: 
	fi: 
	for j from 1 to nops(body) do
		RET:={op(RET), op(1,op(j,body))}:
	od: 
od:
return(RET):
end proc:

#####################################################################

PCAD_ContSet:=proc( pols, var::symbol, $ ) :: set:
local S, i, polys:
description "PCAD_ContSet: Compute the set of contents (with respect to a variable) from a set of polynomials":

if pols::list then 
	polys:=convert(pols,set): 
elif pols::set then 
	polys:=pols:
else 
	error("first input should be a list or set of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
S:={seq( content(op(i,polys),var), i=1..nops(polys))}:
remove(X->X::'constant',S):
end proc:

#####################################################################

PCAD_PrimSet:=proc( pols, var::symbol, $ ) :: set:
local S, i, polys:
description "PCAD_PrimSet: Compute the set of primitive parts (with respect to a variable) from a set of polynomials":

if pols::list then 
	polys:=convert(pols,set): 
elif pols::set then 
	polys:=pols:
else 
	error("first input should be a list or set of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
S:={seq( primpart(op(i,polys),var), i=1..nops(polys))}:
remove(X->X::'constant',S):
end proc:

#####################################################################

PCAD_SetGCD := proc(inpols, $) :: polynom:
local S, Div, i:
description "PCAD_SetGCD: Compute the gcd of a set of polynomials.":

if inpols::list then 
	S:=convert(inpols,set): 
elif inpols::set then 
	S:=inpols:
else 
	error("first input should be a list or set of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
if nops(S)=1 then 
	Div:=op(S): 
else 
	Div:=gcd(op(1,S),op(2,S)): 
fi:
for i from 3 to nops(S) do 
	Div:=gcd(Div,op(i,S)):
od: Div:
end proc:

#####################################################################

PCAD_SetFactors := proc(inpols, $) :: set:
local S, F, i;
description "CADdebug_SetFactors: Compute the factors of a set of polynomials.":

if inpols::list then 
	S:=convert(inpols,set): 
elif inpols::set then 
	S:=inpols:
else 
	error("first input should be a list or set of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
F:=table():
for i from 1 to nops(S) do 
	F[i]:=op(map(X->op(1,X),op(2,factors(op(i,S))))): 
od:
F:={entries(F,'nolist')}:
end proc:

#####################################################################

PCAD_DoIntervalsIntersect := proc( int1, int2, $ )
local L1, L2, R1, R2, innerL, otherR:
description "PCAD_DoIntervalsIntersect: Determines whether two intervals intersect.",
			"Input: Two intervals, each should be a list of two rational numbers.",
			"Output: A boolean determining whether they intersect.":

L1 := op(1, int1): 
R1 := op(2, int1):
L2 := op(1, int2): 
R2 := op(2, int2):
# Determine which is the innermost left end and compare with the others right end.
if L2<L1 then 
	innerL := L1: 
	otherR := R2:
else 
	innerL := L2: 
	otherR := R1: 
fi: 
if otherR < innerL then 
	return(false): 
else 
	return(true): 
fi:
end proc:

#####################################################################

PCAD_RCBtoCube := proc(RCbox::table, R::table) :: list(list):
local Cube, eqs, var:
description "PCAD_RCBtoCube: Convert an object of BoxType from RegularChains to a Cube.":
			"Input: An object of box type (from RegularChains) and a polynomial ring of the same dimension.",  
			"Output: A cube (list of list of rationals) describing the box.":

if R['type']<>'polynomial_ring' then
	error("second argument should be a polynomial ring"):
elif RCbox['type']<>'BoxType' then
	error("first argument should be a box from the RegularChains library"):
fi:
eqs := RCbox['box_bwe'][1]:
Cube :=[]:
for var in ListTools:-Reverse( R['variables'] ) do 
	Cube := [op(Cube), rhs(op(select(X->lhs(X)=var, eqs)))]:
od:
return(Cube):
end proc:

#####################################################################

PCAD_DoCubesIntersect := proc( c1, c2, AllIntersect )
local i, dim, ans, int1, int2:
description "PCAD_DoCubesIntersect: Determines whether two cubes (intervals in many dimensions) intersect.",
				 "The default is that they intersect in ALL dimensions, but optionally can specify just one." ,
			"Input: Two cubes, each should be a list of lists of two rational numbers.",
			"Output: A boolean determining whether the cubes intersect.":

dim := nops(c1):
# Look at each interval in turn
if AllIntersect=true then
	for i from 1 to dim do 
		int1 := c1[i]:
		int2 := c2[i]:
		ans := PCAD_DoIntervalsIntersect( int1, int2 ):
# If we require them all to intersect then any failure may trigger the break and the answer is false.
		if ans=false then 
			break:
		fi:
	od:
else
	for i from 1 to dim do 
		int1 := c1[i]:
		int2 := c2[i]:
		ans := PCAD_DoIntervalsIntersect( int1, int2 ):
# If we only require one to intersect then any success may trigger the break and the answer is true.
		if ans=true then 
			break:
		fi:
	od:
fi:
return(ans):
end proc:

####################################################################

PCAD_IsIntervalInsideAnother := proc( int1, int2, $)
local L1, R1, L2, R2:
description "PCAD_IsIntervalInsideAnother: Determines whether one interval is inside another.",
			"Input: Two intervals, each should be a list of two rational numbers.",
			"Output: A boolean determining whether the first is strictly inside the second.":

L1 := op(1, int1): 
R1 := op(2, int1):
L2 := op(1, int2): 
R2 := op(2, int2):
# Make sure the first is inside the second for true
if L1<=L2 or R1>=R2 then
	return(false):
else
	return(true):
fi:
end proc:

#####################################################################

PCAD_IsCubeInsideAnother := proc( c1, c2, $)
local i, dim, ans, int1, int2:
description "PCAD_IsCubeInsideAnother: Determines whether one cube is inside another.",
			"Input: Two cubes, each should be a list of lists of two rational numbers.",
			"Output: A boolean determining whether the first is strictly inside the second.":

dim := nops(c1):
for i from 1 to dim do 
	int1 := c1[i]:
	int2 := c2[i]:
	ans := PCAD_IsIntervalInsideAnother( int1, int2 ):
	#We require them all to be inside so any failure can trigger false.
	if ans=false then 
		break:
	fi:
od:
return(ans):
end proc:

#####################################################################
###  Section 3: Collins Projection
#####################################################################

PCAD_CCADmat:=proc(S::list, var::symbol, $) :: 'Matrix':
local i, j, k, ell, M;
description "PCAD_CCADmat: Compute the matrix associated with a set of polynomials (see Jirstrand page 15).":

if S=[] then 
	return(Matrix([0])); 
fi:
if {op(map(type,S,polynom))}<>{true} then 
	error("The first input should be a list of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
k:=nops(S): 
ell:=1+max(map(degree,S,var)):
M:=Matrix(k,ell):
for i from 1 to k do 
	for j from 1 to ell do
		M[i,j]:=coeff(S[i],var,ell-j): 
	od: 
od: 
M:
end proc:

#####################################################################

PCAD_CCADpsc:=proc(F::polynom, G::polynom, K::nonnegint, var::symbol, $) :: 'polynom':
local f, g, m, n, i, k, ell, M, count, S;
description "PCAD_CCADpsc: Compute the kth principal subresultant coefficient of two polynomials with respect to a variable.":

if degree(F,var)<degree(G,var) then 
	f:=G: g:=F: 
else 
	f:=F: g:=G: 
fi:
m:=degree(f,var): 
n:=degree(g,var): 
if K>n+1 then 
	error("third input must be an integer less than minimum degree of the input polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
S:=table(): 
count:=1: 
for i from (n-K-1) by -1 to 0 do 
	S[count]:=expand(var^i*f): 
	count:=count+1: 
od:
for i from (m-K-1) by -1 to 0 do 
	S[count]:=expand(var^i*g): 
	count:=count+1: 
od:
S:=[entries(S,'nolist')]: 
M:=PCAD_CCADmat(S,var): 
k:=op(1,[LinearAlgebra:-Dimensions(M)]); 
ell:=op(2,[LinearAlgebra:-Dimensions(M)]): 
if k=ell then 
	return(LinearAlgebra:-Determinant(M)): 
else 
	return( LinearAlgebra:-Determinant( LinearAlgebra:-DeleteColumn(M,[op({seq(k..ell)} minus {ell-K})]) )):
fi: 
end proc:

#####################################################################

PCAD_CCADpscset:=proc(f::polynom, g::polynom, var::symbol, $) :: 'set':
local n, j:
description "PCAD_CCADpscset: Compute the set of principal subresultant coefficients of two polynomials with respect to a variable.":

n:=min(degree(f,var),degree(g,var)): 
{seq( PCAD_CCADpsc(f,g,j,var), j=0..n)} minus {0}:
end proc:

#####################################################################

PCAD_CCADreductum:=proc(poly::polynom, var::symbol, deg::nonnegint, $) :: 'polynom':
description "PCAD_CCADreductum: Compute the reductum of a polynomial.":

if deg>degree(poly,var) then 
	error("reductum degree is bigger than polynomial degree. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
convert(series( poly,var,deg+1 ),polynom):
end proc:

#####################################################################

PCAD_CCADreductumset:=proc(poly::polynom, var::symbol, lvars::set, $) :: 'set':
local i, j, deg, rlist, cflist, rset, count, fin, tvar, res, crw:
description "PCAD_CCADreductumset: Compute the reductum set of a polynomial.  Note that this is the simplified reductum set as required for Collins projection.":

[var,op(lvars)]:
deg:=degree(poly,var): 
cflist:=ListTools:-Reverse( [seq( coeff(poly,var,i), i=0..deg)] ):
rlist:=ListTools:-Reverse( [seq( PCAD_CCADreductum(poly,var,j), j=0..deg)] ):
crw:=op(1,cflist):
fin:=false:
if crw::'constant' or nops(lvars)=1 then 
	rset:=rlist[1..1]: 
	fin:=true:
fi:
count:=2:
while fin=false do 
	if count>=nops(cflist) then 
		fin:=true: 
		rset:=rlist: 
	else 
		tvar:=indets(crw,'name')[1]:
		res:=resultant(crw, op(count,cflist), tvar): 
		if res=0 then 
			crw:=crw:
			count:=count+1:
		elif res::'constant' then 
			rset:=[op(1..count-1,rlist)]: 
			fin:=true:
		else 
			crw:=res:
			count:=count+1:
		fi:
	fi:
od:
rset:={op(rset)} minus {0}:
end proc:

#####################################################################

PCAD_CCADProj1:=proc( polyset::set, var::symbol, lvars::list, $ ) :: 'set':
local i, entry, redset, S, count:
description "PCAD_CCADProj1: Compute the first part of the Collins projection operator.":

if map(X->type(X,polynom),polyset)<>{true} then 
	error("input should be set of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"); 
fi:
S:=table(): 
count:=1:
for i from 1 to nops(polyset) do 
	redset:=PCAD_CCADreductumset(op(i,polyset), var, convert(lvars,'set')):
	for entry in redset do
		S[count]:=PCAD_LdCf(entry,var): count:=count+1:
		S[count]:=op(PCAD_CCADpscset(entry,diff(entry,var),var)): count:=count+1:
	od:
od:
S:={entries(S,'nolist')}:
S:=remove(X->X::'constant',S):
S:=PCAD_RemoveConstantMultiples(S):
end proc:

#####################################################################

PCAD_CCADProj2:=proc( polyset::set, var::symbol, lvars::list, $ ) :: 'set':
local i, j, entry, ent, redset1, redset2, S, count:
description "PCAD_CCADProj2: The second part of the Collins projection operator.":

if map(X->type(X,polynom),polyset)<>{true} then 
	error("input should be set of polynomials. This is an unexpected error message - please report to M.England@bath.ac.uk"); 
fi:
if nops(polyset)=1 then 
	return({}): 
fi: 
S:=table(): count:=1:
for i from 1 to nops(polyset) do
	redset1:=PCAD_CCADreductumset(op(i,polyset),var, convert(lvars,'set') ): 
	for j from i+1 to nops(polyset) do
		redset2:=PCAD_CCADreductumset(op(j,polyset),var, convert(lvars,'set')):
		for entry in redset1 do 
			for ent in redset2 do
				S[count]:=op(PCAD_CCADpscset(entry,ent,var)): 
				count:=count+1:
			od: 
		od:
	od:
od:
S:={entries(S,'nolist')}:
S:=remove(X->X::'constant',S):
S:=PCAD_RemoveConstantMultiples(S):
end proc:

#####################################################################

PCAD_CCADProj:=proc( polyset::set, var::symbol, lvars::list, {SF::boolean:=true}, $ ) :: 'set':
local pset, S1, S2, S:
description "PCAD_CCADProj: The full Collins projection operator.":

userinfo(3,{'PCAD_CCADProj', 'CADProjection'},"now using Collins projection with respect to", var):
if SF=true then 
	pset:=PCAD_SFBasis(polyset,var): 
	userinfo(3,{'PCAD_CCADProj', 'CADProjection'},"square free basis consisting of ", nops(pset), "polynomials"):
else 
	pset:=polyset: 
	userinfo(3,{'PCAD_CCADProj', 'CADProjection'},"set of ", nops(pset), "polynomials"):
fi:
userinfo(4,{'PCAD_CCADProj', 'CADProjection'},"the polynomials are", pset):
S1:=PCAD_CCADProj1(pset,var,lvars):
userinfo(3,{'PCAD_CCADProj', 'CADProjection'},"first operator added ", nops(S1)):
S2:=PCAD_CCADProj2(pset,var,lvars):
userinfo(3,{'PCAD_CCADProj', 'CADProjection'},"second operator added ", nops(S2)):
S:=S1 union S2:
S:=PCAD_RemoveConstantMultiples(S):
end proc:

#####################################################################

PCAD_CCADProjPolys:=proc( polyset::set, vars::list, {SF::boolean:=true}, $ ) :: 'set':
local i, pset, ret, cont, dim:
description "PCAD_CCADProjPolys: Compute the set of projection polynomials for a given set of polynomials using the Collins projection operator and with respect to a list of variables in descending order.":

dim:=nops(vars):
pset:=table(): 
pset[0] := PCAD_PrimSet(polyset,vars[1]):
pset[0] := PCAD_SFBasis(pset[0], vars[1]):
pset[0] := PCAD_SetFactors(pset[0]):
cont:=PCAD_ContSet(polyset,vars[1]): 
for i from 1 to dim-1 do
	pset[i]:=PCAD_CCADProj( pset[i-1], vars[i], vars[i+1..dim] ) union cont:
	cont:=PCAD_ContSet(pset[i],vars[i+1]):
	pset[i]:=PCAD_PrimSet(pset[i], vars[i+1]):
	if i=nops(vars)-1 then 
		if SF=true then 
			pset[i]:=PCAD_SFBasis(pset[i],vars[i+1]): 
		fi:
	else 
		pset[i]:=PCAD_SFBasis(pset[i],vars[i+1]):
	fi:
	pset[i]:=PCAD_SetFactors(pset[i]):
od:
ret:=map(X->op(X),[entries(pset,'nolist')]):
ret:=PCAD_RemoveConstantMultiples(ret):
end proc:

#####################################################################
### Section 4: McCallum Projection
#####################################################################

PCAD_McCallumProj:=proc( polys::set, var::symbol, lvars::list, $) :: set:
local Polys, cont, Pol, clist, Pset1, Pset2, Pset, i, j:
description "PCAD_McCallumProj: The McCallum projection operator.":

Polys:=PCAD_PrimSet(polys,var): 
cont:=PCAD_ContSet(polys,var);
Polys:=PCAD_SFBasis(Polys,var): 
userinfo(3,{'PCAD_McCallumProj', 'CADProjection'},"now using McCallum projection with respect to", var):
userinfo(3,{'PCAD_McCallumProj', 'CADProjection'},"square free basis of primitive part consisting of ", nops(Polys), "polynomials"):
userinfo(4,{'PCAD_McCallumProj', 'CADProjection'},"the polynomials are", Polys):
userinfo(4,{'PCAD_McCallumProj', 'CADProjection'},"content ", cont, " was added to output."):
Pset1:=table():
for i from 1 to nops(Polys) do 
	Pol:=op(i,Polys): 
	clist:=PCAD_McCcoeffset(Pol,var, convert(lvars,'set')):
	Pset1[i]:=discrim(Pol,var),op(clist):
	userinfo(4,{'PCAD_McCallumProj', 'CADProjection'},"added discriminant of polynomial", i, "and first", nops(clist), "coeffs"): 
od:
Pset2:=table():
for i from 1 to nops(Polys) do 
	for j from i+1 to nops(Polys) do 
		Pset2[i,j]:=resultant( op(i,Polys), op(j,Polys), var):
		userinfo(4,{'PCAD_McCallumProj', 'CADProjection'},"added resultant of polynomials ", i, "and ", j):
	od:
od:
Pset:={op(cont), entries(Pset1,'nolist'), entries(Pset2,'nolist')}:
Pset:=remove(X->type(X,constant),Pset): 
userinfo(3,{'PCAD_McCallumProj', 'CADProjection'}, nops(Pset), "non constant projection polynomials found"):
Pset:=PCAD_RemoveConstantMultiples(Pset): 
userinfo(3,{'PCAD_McCallumProj', 'CADProjection'},"removed constants and constant multiples to leave", nops(Pset)):
Pset:
end proc:

#####################################################################

PCAD_McCcoeffset:=proc( poly::polynom, var::symbol, lvars::set, $ ) :: list:
local i, crw, deg, cflist, tvar, res:
description "PCAD_McCcoeffset: Compute the set of coefficients required by the McCallum projection operator.":

deg:=degree(poly,var): 
cflist:=ListTools:-Reverse( [seq( coeff(poly,var,i), i=0..deg)] ):
crw:=op(1,cflist): 
if crw::'constant' then 
	return([]): 
elif nops(lvars)=1 then
	return([crw]):
fi:
for i from 2 to nops(cflist) do 
	tvar:=indets(crw,'name')[1]:
	res:=resultant(crw,op(i,cflist), tvar):
	if res=0 then 
		crw:=crw:
	elif res::'constant' then 
		return([op(1..i,cflist)]):
	else 
		crw:=res:
	fi:
od:
cflist:
end proc:

#####################################################################

PCAD_McCProjPolys:=proc( polyset::set, vars::list, {SF::boolean:=true}, $ )::'set':
local i, pset, ret, cont, dim:
description "PCAD_McCProjPolys: Compute the set of projection polynomials for a given set of polynomials using the McCallum projection operator and with respect to a list of variables in descending order.":

dim:=nops(vars):
pset:=table(): 
pset[0] := PCAD_PrimSet(polyset,vars[1]):
pset[0] := PCAD_SFBasis(pset[0], vars[1]):
pset[0] := PCAD_SetFactors(pset[0]):
cont:=PCAD_ContSet(polyset,vars[1]): 
for i from 1 to dim-1 do
	pset[i]:=PCAD_McCallumProj( pset[i-1], vars[i], vars[i+1..dim]  ) union cont:
	cont:=PCAD_ContSet(pset[i],vars[i+1]):
	pset[i]:=PCAD_PrimSet(pset[i], vars[i+1]):
	if i=nops(vars)-1 then 
		if SF=true then 
			pset[i]:=PCAD_SFBasis(pset[i],vars[i+1]): 
		fi:
	else 
		pset[i]:=PCAD_SFBasis(pset[i],vars[i+1]):
	fi:
	pset[i]:=PCAD_SetFactors(pset[i]):
od:
ret:=map(X->op(X),[entries(pset,'nolist')]):
ret:=PCAD_RemoveConstantMultiples(ret):
end proc:

#####################################################################
### Section 5: Lifting
#####################################################################

PCAD_IsCellZeroDim:=proc( cellindex::list, $ ) :: boolean:
local alleven:
description "PCAD_IsCellZeroDim: Check to see if a cell is zero dimensional by examining the index.":

alleven:=convert(map(X->is(X::'even'),cellindex),set): 
if alleven={true} then 
	return(true): 
elif alleven={false} or alleven={false,true} then 
	return(false): 
else 
	error("could not determine parity. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
end proc: 

#####################################################################

PCAD_IsCellFullDim:=proc( cellindex::list, $ ) :: boolean:
local allodd:
description "PCAD_IsCellFullDim: Check to see if a cell is of full dimension by examining the index.":

allodd:=convert(map(X->is(X::'odd'),cellindex),set): 
if allodd={true} then 
	return(true): 
elif allodd={false} or allodd={false,true} then 
	return(false): 
else 
	error("could not determine parity. This is an unexpected error message - please report to M.England@bath.ac.uk"):
fi:
end proc: 

#####################################################################

PCAD_MinimalDelineatingPolynomial:=proc(f::polynom,alpha::list, R, $) :: polynom:
local k, j, vars, mdeg, t, Dp, S, vchoice:
description "PCAD_MinimalDelineatingPolynomial: Find the minimal delineating polynomial following Brown's 2005 technical note.":

vars:=R['variables']: 
if RegularChains:-MainVariable(f,R)<>vars[1] then 
	error("incorrect polynomial ring in input. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
mdeg:=degree(f,vars): 
t:=0: 
k:=1:
while t=0 do 
	if k>mdeg then 
		error("the minimal delineating polynomial was not found. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
	fi:
	vchoice:=map(X -> seq(X+0*j,j=1..k),vars):
	Dp:=combinat:-choose(vchoice,k):
	Dp:=map(X->diff(f,X), Dp):
	S:=remove(X->is(subs(alpha,X)=0)=true, Dp):
	if S<>[] then 
		t:=k: 
	fi:
	k:=k+1:
od: 
PCAD_SetGCD(S):
end proc:

#####################################################################

PCAD_CanRCHaveSolutionAtSP := proc( rc, SPrc, C::list(list(rational)), R, $) :: boolean:
local vars, eqs, speqs, ineqs, lR, lrc, rootboxes, rootcubes, TheCube, ans, gap, Poss, def:
description "PCAD_CanRCHaveSolutionAtSP: Given a regular chain and a sample point of one dimension lower (regular chain and cube) determine if it can have a solution.",
			"Input: A regular chain, a sample point stored as regular chain and cube (list of list of rationals) and a polynomial ring.  The ring and chain should have the same dimension and the cube one dimension lower.",  
			"Output: A boolean indicating whether a solutions is possible.  (Note: we are only checking the lower dim part of zero dim chains so no solution in the top dimension is guarenteed).":


# Extracting lower dimension information from system.
vars := R['variables'][2..-1]:
eqs := RegularChains:-Equations(rc, R)[2..-1]:
eqs :=  map(X -> PCAD_MakeMonic(X, R), eqs):
ineqs := RegularChains:-Inequations(rc, R):
ineqs := remove(has, ineqs, R['variables'][1]):
# Quick check for easy exit - if the lower dim parts are the same then return true.
speqs := map(X -> PCAD_MakeMonic(X, R), RegularChains:-Equations(SPrc, R)):
if convert(eqs, 'set') = convert(speqs, 'set') then
	return(true):
fi:
# Forming lower dimension rc
lR := RegularChains:-PolynomialRing( vars ):
lrc := RegularChains:-Triangularize( eqs, ineqs, lR ):
# Output of Triangularize should be a list with one RC given where eqs came from
if nops(lrc)=1 and lrc[1]['type']='regular_chain' then
	 lrc := op(lrc):
else
	error("Output of Triangularize not as expected.  Please report to M.England@bath.ac.uk"):
fi:
# If chain is not zero dimensional then do not check, just return true.
if RegularChains:-ChainTools:-IsZeroDimensional(lrc, lR)<>true then
	return(true):
fi:
# We look at the roots of the system and keep refining until none of the cubes intersect the BB or just one is inside it. 
ans := FAIL:
gap := 1/10000:
while ans=FAIL do
	# Go for smaller bounding boxes.
	gap := gap*(1/100):
	# Find a list of boxes describing the real roots for the lower dimension system.
	rootboxes := RegularChains:-SemiAlgebraicSetTools:-RealRootIsolate( eqs, [], [], ineqs, lR, 'abserr'=gap):
	# Converting regular chain box type to basic cube type used in ProjectionCAD. 
	rootcubes := map(X -> PCAD_RCBtoCube(X, lR), rootboxes):
	# List of those cubes which intersect the BB (in all dimensions). 
	Poss := select(X -> PCAD_DoCubesIntersect(X, C, true)=true, rootcubes):
	if Poss=[] then
		# If none of the cubes intersect then this system has no compatible root, so return false.
		return(false):
	elif nops(Poss)<>1 then
		# If more than one of the cubes is intersecting the BB then we need to refine more and go again.
		next:
	else
		TheCube := op(Poss):
		def := PCAD_IsCubeInsideAnother(TheCube, C):
		if def=true then
			# The system has a single compatible root so return true.
			return(true):
		else
			# There is one possible root but we are not certain - refine more.
			next:
		fi:
	fi:
od:
end proc:

#####################################################################

PCAD_CanRSHaveSolutionAtSP := proc( rs, SPrc, C::list(list(rational)), R, $) :: boolean:
local rc, ans:
description "PCAD_CanRSHaveSolutionAtSP: Given a regular system and sample point of one dim lower (regular chain and cube) determine if it can have a solution in the cube.",
			"Input: A regular system, a sample point (regular chain and cube (list of list of rationals)) and a polynomial ring.  The ring and system should have the same dimension and the sample point one dimension lower.",  
			"Output: A boolean indicating whether a solutions is possible.  (Note: we are only checking the lower dim part so no solution in the top dimension is guarenteed).":

# We simple extract the RC and see if that is compatible.  
	# We are ignoring the inequations.  This means we may return true when the answer is really false (but the opposite is not possible).
	# Hence the worst case is including an unrequired poly in lifting set and thus having more cells than required.  Output is still correct though.
rc := RegularChains:-ConstructibleSetTools:-RepresentingChain(rs, R):
ans := PCAD_CanRCHaveSolutionAtSP(rc, SPrc, C, R):
return(ans):
end proc:

#####################################################################

PCAD_MakeCoprimeOverPoint := proc( RPrc::table, SPrc::table, SPcube::list, pset::list, R, $ ) :: 'list':
local i, j, RC, out, ret, new, RSs, rs, tmprc, mvar, eqs, inits, Hinits, RCeqs:
description "PCAD_MakeCoprimeOverPoint: Given a real point and a set of polynomials, produce another set of polynomials with equivalent zeros which are coprime over the real point.",
			"Input: The sample point encoded as a regular chain and cube (list of lists of pairs of rationals) isolating a single root and the RC defining the real point part.",
				"Also, a polynomial ring is provided of the same dimension and a set of polynomials (with same main variable).",  
			"Output: A list of polynomials.":

# Trivial case first (shouldn't actually occur I think).
if pset=[] then 
	return(pset): 
fi:
mvar := R['variables'][1]:
RC := RPrc:
RCeqs := convert(map(X -> PCAD_MakeMonic(X,R), RegularChains:-Equations(RC, R)), 'set'): 
ret := []: 
for i from 1 to nops(pset) do
	#We repeatedly look for the solutions to the equations, such that (a) the base rc is satisfied and (b) the equations in ret (from previous steps) are not zero. 
	out := RegularChains:-Triangularize( [op(i, pset)], ret, RC, R ):
	# There are two types of output from Triangularize
	# First consider list of regular chains.  
	if out=[] then
		# If an empty list then Triangulaize has identified that no solutions are possible and thus returns an empty list.
		ret := ret:
	elif out::list then 
		if convert(map(X -> X['type'], out), 'set')<>{'regular_chain'} then
			error("output from Triangularize was not as expected. This is an unexpected error message - please report to M.England@bath.ac.uk"):
		fi:
		# This next line removes any regular chains whose first equation is not in the main variable.  
			#This could occur for example if the lower point nullifies the higher equation so there is a solution for all mvar. 
		out := select(X -> RegularChains:-MainVariable( RegularChains:-Equations(X, R)[1], R)=mvar, out):
		# Now we check if the rc produced with Triangularize is compatible in lower dimension with the cube.
		# Note: Previosuly I only checked when more than one (a split) but examples show even with one if may be incompatible (in which case the polynomial is not required).
		out := select(X -> PCAD_CanRCHaveSolutionAtSP(X, SPrc, SPcube, R)=true, out):
		if out=[] then
			ret := ret:
		else
			# Extract polynomial in main variable (if the correct level).
			for j from 1 to nops(out) do
				tmprc := op(j,out): 
				new:=expand(RegularChains:-Equations(tmprc,R)[1]): 
				if RegularChains:-MainVariable(new,R)=mvar then
					ret := [op(ret),new]:
				fi:
			od:
		fi:
	# Second consider when output is a constructable set.
	elif out::table and out[type]='constructible_set' then
		RSs:=out['list_regular_system']: Display(%,R);
		if RSs=[] then
			# If an empty list then Triangulaize has identified that no solutions are possible and thus returns an empty list.
			ret := ret:
		else
			if convert(map(X -> X['type'], RSs), 'set')<>{'regular_system'} then
				error("output from Triangularize was not as expected. This is an unexpected error message - please report to M.England@bath.ac.uk"):
			fi:
			# This next line removes any regular systems whose first equation is not in the main variable.  
				#This could occur for example if the lower point nullifies the higher equation so there is a solution for all mvar. 
			RSs := select(X -> RegularChains:-MainVariable( RegularChains:-Equations(RegularChains:-ConstructibleSetTools:-RepresentingChain(X, R), R)[1], R)=mvar, RSs):
			# Now we check if the RSs are compatible in lower dimension with the cube.
				# Note: Previosuly I only checked when more than one (a split) but examples show even with one if may be incompatible (in which case the polynomial is not required).
			RSs := select(X -> PCAD_CanRSHaveSolutionAtSP(X, SPrc, SPcube, R)=true, RSs):
			if RSs=[] then
				ret := ret:
			else
				# Extract polynomial in main variable and keep.
				for j from 1 to nops(RSs) do
					rs := op(j,RSs): 
					eqs := map(X -> PCAD_MakeMonic(X,R), RegularChains:-Equations(rs['regular_system_chain'], R)):
					inits := convert(map(X -> PCAD_MakeMonic(RegularChains:-Initial(X,R),R), ret), 'set'):
					Hinits := select(X -> has(inits, X)=true, eqs):
					Hinits := select(X -> has(RCeqs, X)=false, Hinits):
					if Hinits=[] then
						new := eqs[1]: 
						if RegularChains:-MainVariable(new,R)=mvar then
							ret := [op(ret),new]:
						fi:
					# else this is the special case where the initial of one of the proj poly is forced zero, but we are not in the cell where that happens.  Hence we ignore this case. 
					fi:
				od:
			fi:
		fi:
	else 
		error("output from Triangularize was not as expected. This is an unexpected error message - please report to M.England@bath.ac.uk"):
	fi:
od:
return( convert(convert(map(X->PCAD_MakeInteger(X),ret),'set'),'list') ):
end proc:

#####################################################################

PCAD_WhichRCInSplitHasSolutionAtSP := proc( LRCs::list, Cube::list, R, $) :: posint:
local i, dim, lowerRCeqs, lowerRC, lowerR, lowerCube, MainPolys, out, IndexList, tRoots, answer, OMswitch:
description "PCAD_WhichRCInSplitHasSolutionAtSP: Given a list of regular chains and an isolating cube, identify which entry in the list has a zero in the cube, returning that index.",
			"Input: A list of regularchains, a cube (list of lists of pairs of rationals) and a polynomial ring.  The cube is one dimension lower than the others.",
			"Output: A positive integer denoting which of the regular chains has a solutions within the cube.",
				"Assumption:  There is an assumption here that all the rcs came from by splitting a single one, and thus that the cube isolates a single root.":

# Basic error checking.
if {op(map(X->X[type],LRCs))}<>{'regular_chain'} then 
	error("first input should be a list of regular chains. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
if R[type]<>'polynomial_ring' then 
	error("third input should be a polynomial ring. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
# We use internal regular chains commands so need opaque modules off.
if kernelopts('opaquemodules')=true then 
	OMswitch := 0:
	kernelopts('opaquemodules'=false): 
fi:
# Identify lower dim parts.
lowerRCeqs := map(X -> PCAD_SortEqsRC(RegularChains:-Equations(X,R),R)[1..-2], LRCs):
lowerRCeqs := map(X -> map(Y -> PCAD_MakeMonic(Y,R), X), lowerRCeqs):
if nops(convert(lowerRCeqs,set))<>1 then
	error("the regular chains inputted differ at a lower level. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
else 
	lowerRCeqs := op(1,lowerRCeqs): 
fi:
dim := nops(R['variables']):
# Isolate zeros to see if they match.
if dim=1 then 
	lowerRC := RegularChains:-ChainTools:-Empty(R):
	lowerCube := []:
	MainPolys := map(X->RegularChains:-Equations(X,R)[1], LRCs): 
	out := RegularChains:-TRDisolate_real_zeros( [lowerRC,lowerCube], MainPolys, R ):
else
	lowerR := RegularChains:-PolynomialRing( R['variables'][2..-1] ): 
	lowerRC := RegularChains:-ChainTools:-Chain( lowerRCeqs, RegularChains:-ChainTools:-Empty(lowerR), lowerR):
	lowerCube := Cube[1..-2]:
	MainPolys := map(X->RegularChains:-Equations(X,R)[1], LRCs): 
	out := RegularChains:-TRDisolate_real_zeros( [lowerRC,lowerCube], MainPolys, R ):
fi:
IndexList := op(2,out): 
tRoots := map(X->op(-1,X),op(1,out)): 
tRoots := map(X->PCAD_DoIntervalsIntersect(X,Cube[-1]), tRoots):
# Should only be one answer.
if nops(select(X->X=true,tRoots))<>1 then 
	error("the isolating cubes overlap. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
# Identify answer.
for i from 1 to nops(tRoots) do 
	if op(i,tRoots)=true then 
		answer := IndexList[i]: 
	fi: 
od:
# Leave user settings as we find them.
if OMswitch=0 then 
	kernelopts('opaquemodules'=true): 
fi:
return(answer):
end proc:

#####################################################################

PCAD_MakeSquareFreeOverPoint := proc( rc::table, extraeqs::list, SPcube::list, pset::list, R, $ ) :: list:
local i, j, ret, out, newpol, possRCs, rightindex, newrc, fin, pol, tail, Check: 
description "PCAD_MakeSquareFreeOverPoint: Given a real point encoded by regular chain and isolating cube, and a set of polynomials, produce a set of polynomials with equivalent zeros which are squarefree over the real point.",
			"Input: The real point is encoded as a regular chain and cube (list of lists of pairs of rationals) isolating a single root.  A polynomial ring is provided of the same dimension and a set of polynomials (with same main variable).",
				"The extraeqs input is outdated and not used any more.",
			"Output: A list of polynomials, squarefree over the point.":

# Trivial case first.
if pset=[] then 
	return(pset): 
fi:
ret := []:
for i from 1 to nops(pset) do
	if RegularChains:-MainVariable(pset[i],R)<>R['variables'][1] then 
		ret := [op(ret),pset[i]]:
	else
		pol := pset[i]:
		fin := false:
		while fin=false do
			# To use RegularChains:-ChainTools:-SquarefreeFactorization on a polynomial it must be regular so we check this.
			Check := RegularChains:-IsRegular(RegularChains:-Initial(pol, R), rc, R):
			if Check=true then
				out := RegularChains:-ChainTools:-SquarefreeFactorization( pol, R['variables'][1], rc, R, 'method'='src'):
				fin := true:
			else
				# If not regular then consider tail (as leading coefficient will vanish).
				tail := RegularChains:-Tail(pol, R):
				if tail::'constant' then
					#If we have lost the main variable then we do not need this polynomial.
					out := false:
					fin := true:
				elif RegularChains:-MainVariable(pol,R)=RegularChains:-MainVariable(tail,R) then
					# Run again on tail
					pol := tail:
				else
					#If we have lost the main variable then we do not need this polynomial.
					out := false:
					fin := true:
				fi:
			fi:
		od:
		if out=false then
			# We didn't need this polynomial
			ret := ret:
		else
			if nops(out)<>1 then
				# This means the regular chain has been split, but we actually only care about one of the roots and thus one rc in the split so we identify which.
				possRCs := map(X -> X[-1],out): 
				possRCs := map(X -> RegularChains:-ChainTools:-Chain( PCAD_SortEqsRC( [op(extraeqs),op(RegularChains:-Equations(X,R))] ,R), RegularChains:-ChainTools:-Empty(R), R ), possRCs):
				rightindex := PCAD_WhichRCInSplitHasSolutionAtSP( possRCs, SPcube, RegularChains:-PolynomialRing( R['variables'][2..-1] ) ): 
				out := out[rightindex]:
			else 
				out := op(out):
			fi:
			newrc := op(2,out): 
			if RegularChains:-ChainTools:-IsIncluded(newrc,rc,R)<>true then 
				error("SquarefreeFactorization output is not as expected. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
			fi:
			out:=op(1,out):
			newpol:=1:
			for j from 1 to nops(out) do
				newpol := newpol*op(1,op(j,out)):
			od: 
			ret := [op(ret), expand(newpol) ]: 
		fi:
	fi:
od:
return(ret):
end proc: 

#####################################################################

PCAD_SPtoRootOf:=proc(rc::table, cube::list, R, $) :: list:
local i, dim, eqs, vars, pt:
description "PCAD_SPtoRootOf: Convert a sample point encoded by regular chain and isolating cube into a list of substitutions giving the variables as algebraic numbers.",
			"Input: A regular chain and cube (list of lists of pairs of rationals) which together encode a real point.  Also a polynomial ring of the same dimension.",
			"Output: A list of equations setting the variables to RootOf expressions.":

if rc[type]<>'regular_chain' then 
	error("first input should be regular chain. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
if R[type]<>'polynomial_ring' then 
	error("third input should be a polynomial ring. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi:
dim:=nops(cube): 
vars:=[op(1..dim, ListTools:-Reverse( R['variables'] ))]:
eqs:=ListTools:-Reverse( RegularChains:-Equations(rc,R) ): 
if nops(eqs)<>nops(cube) then 
	error("cannot isolate a solution of regular chain. This is an unexpected error message - please report to M.England@bath.ac.uk"): 
fi: 
pt:=[]:
for i from 1 to dim do 
	if op(i,cube)[1]=op(i,cube)[2] then 
		pt:=[op(pt), vars[i]=op(i,cube)[1]]:
	else 
		pt:=[op(pt), vars[i]=RootOf( subs(pt,eqs[i]), op(i,cube)[1]..op(i,cube)[2]) ]:
	fi:
od:
end proc:

#####################################################################

PCAD_SortCadCell:=proc( c1, c2, $ ) :: boolean:
local tmp:
description "PCAD_SortCadCell: A rule to decide the order of cadcells according to their index. For use as a rule in sort.":

tmp:=[op(1,c1),op(1,c2)]:
if sort( tmp ) = tmp then 
	return(true):
else 
	return(false):
fi:
end proc:

#####################################################################

PCAD_GenerateStack := proc( cell::list, in_F::list, R, output, $ ) :: 'list':
local i, F, OMswitch, cellIndex, cellSP, cellSPrc, cellSPrcEqs, cellSPcube, RPrc, RPrcEqs, retstack: 
description "PCAD_GenerateStack: Generate a stack over a cell with respect to a given polynomials. This uses the internal RegularChains command but first preconditions the polynomials so they are of the right format.",
			"Input: A cell in list or listwithrep format, a list of polynomials and a polynomial ring.  The ring and polynomials have mvar one higher than the cell.  Also, the output format is specified.",  
			"Output: The stack encoded as a list of cells in the states output format.":

# We will need internal RegularChains commands so need opaquemodules off if it isn't already.
if kernelopts('opaquemodules')=true then 
	OMswitch := 0: 
	kernelopts('opaquemodules'=false): 
fi:
# Extract info from cell
F := expand(in_F):
cellIndex:=cell[1]:
cellSP := cell[-1]:
cellSPrc := cellSP[1]:
cellSPrcEqs := ListTools:-Reverse( RegularChains:-Equations(cellSPrc,R) ): 
cellSPcube := cellSP[2]:
# We need to ensure the polynomials in F are coprime and squarefree over the cell.  
	# We have commands to ensure they are coprime and squarefree over a real point.
	# We enforce these over the point in the subspace of R^n implied by the cell - (i.e. the zero dim part if any).
	# Given delineability from the ProjectionCAD theory this will be sufficient.  
RPrcEqs:=[]:
for i from 1 to nops(cellIndex) do 
	if cellIndex[i]::'even' then
		RPrcEqs := [ op(RPrcEqs), cellSPrcEqs[i] ]:
	fi:
od:
####extra := convert( convert(cellSPrcEqs, 'set') minus convert(RPrcEqs, 'set'), 'list'):
if RPrcEqs=[] then 
	if not PCAD_IsCellFullDim(cellIndex) then  
		error("the regular chain was not as expected. This is an unexpected error message - please report to M.England@bath.ac.uk"):
	fi:
else 
	RPrc := RegularChains:-ChainTools:-Chain( PCAD_SortEqsRC(RPrcEqs,R), RegularChains:-ChainTools:-Empty(R), R): 
	F := PCAD_MakeCoprimeOverPoint( RPrc, cellSPrc, cellSPcube, F, R ): 
	F := PCAD_MakeSquareFreeOverPoint( RPrc, [], cellSPcube, F, R ):
fi:
# Use set of factors (without this it seems a bug in interval arithmetic of regular chains can be triggered).  
F := convert( PCAD_SetFactors(F), 'list'):
# Generate stack using RegularChains commands.
if output='list' then
	retstack := RegularChains:-TRDgenerate_stack( [ cell, cell[1] ], F, R ):
else
	retstack := RegularChains:-TRDgenerate_stack_general( [ cell, cell[1] ], F, R):
fi:
# Leave user options as we found them.
if OMswitch=0 then 
	kernelopts('opaquemodules'=true): 
fi:
# Return the stack.
retstack:
end proc:

#####################################################################

PCAD_ProjCADLift:=proc( ProjPolys::set, vars::list, method, finalCAD, retcad, out, failure, $ )
local i, j, k, num, OMswitch, dim, R, contset, pset, cad, tmpcad, tmppset, cell, SPrc, SPcube, cellIndex, alpha, falpha, DP, ICZD:
description "PCAD_ProjCADLift: Given a list of projection polynomials and a list of variables, compute a CAD.":

if kernelopts('opaquemodules')=true then 
	OMswitch:=0: 
	kernelopts('opaquemodules'=false): 
fi:
dim:=nops(vars): 
contset:={}:
for i from 1 to dim do 
	R[i]:=RegularChains:-PolynomialRing( vars[-i..-1] ):
	num:=dim+1-i:
	pset[num]:=select(has, remove(has,ProjPolys,vars[1..i-1]), vars[i]): 
	pset[num]:=pset[num] union contset:
	contset:=PCAD_ContSet(pset[num], vars[i]): 
	pset[num]:=PCAD_PrimSet(pset[num], vars[i]): 
	pset[num]:=convert(pset[num],'list'): 
	pset[num]:=PCAD_SetFactors(pset[num]): 
	pset[num]:=convert(pset[num],'list'): 
	pset[num]:=sort(pset[num],'length'):
	if pset[i]=[] then WARNING("there are no projection polynomials in %1", R[i]['variables']): fi:
od:
pset[1]:=remove(X -> sturm(X, vars[-1], -infinity, infinity) = 0, pset[1]):
if out='list' then
	cad[1]:=RegularChains:-TRDgenerate_stack([], pset[1], R[1]):
else
	cad[1]:=RegularChains:-TRDgenerate_stack_general( [], pset[1], R[1] ):
fi:
cad[1]:=map(X->X[1], cad[1]):
userinfo(2, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"produced CAD of", R[1]['variables'], "-space with", nops(cad[1]), "cells"):
if retcad=1 then
	return(cad[1]):
fi:
for i from 2 to dim do 
	tmpcad:=table():
	for j from 1 to nops(cad[i-1]) do 
		cell:=cad[i-1][j]: 
		cellIndex:=cell[1]: 
		SPrc:=op(1,op(-1,cell)): 
		SPcube:=op(2,op(-1,cell)): 
		alpha:=PCAD_SPtoRootOf( SPrc,SPcube,R[i-1] ):
		tmppset:=[]:
		for k from 1 to nops(pset[i]) do
			falpha:=subs(alpha,op(k,pset[i])): 
			falpha:=is(falpha=0): 
			if falpha=true then
				userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"a projection polynomial was nullified on cell", cellIndex):
				ICZD:=PCAD_IsCellZeroDim( cellIndex ):
				if i=dim then
					if finalCAD='SI' then 
						tmppset:=[op(tmppset)]:
						userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"the nullified polynomial was discarded since this is the final lift and only sign invariance is required."): 
					else 
						if ICZD=true then 
							DP:=PCAD_MinimalDelineatingPolynomial(op(k,pset[i]),alpha,R[i]):
							if DP::constant then
								tmppset:=[op(tmppset)]:
								userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"the nullified polynomial was discarded since cell", cellIndex, "is zero-dim and the minimal delineating polynomial is a constant."):
							else
								tmppset:=[op(tmppset), DP]:
								userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "cell", cellIndex, "is zero-dim so a delineating polynomial was used."):
								userinfo(4, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "the delineating polynomial used is", DP):
							fi:
						else 
							userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "there is nullification on cell", cellIndex, "which has dim>0. Hence the CAD cannot be guaranteed order-invariant."):
							userinfo(4, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "the polynomial", op(k,pset[i]), "is nullified by", alpha, "which is a cell of dim>0."):
							if failure='err' then
								error("The input is not well-oriented (there is nullification on cell %1).  The output cannot be guaranteed correct.", cellIndex):
							elif failure='giveFAIL' then
								return(FAIL):
							else
								WARNING("The input is not well-oriented (there is nullification on cell %1).  The output cannot be guaranteed correct.", cellIndex):
								tmppset:=[op(tmppset)]:
							fi:
						fi:
					fi:
				else
					if method='McCallum' then
						if ICZD=true then 
							DP:=PCAD_MinimalDelineatingPolynomial(op(k,pset[i]),alpha,R[i]):
							if DP::constant then
								tmppset:=[op(tmppset)]:
								userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"the nullified polynomial was discarded since cell", cellIndex, "is zero-dim and the minimal delineating polynomial is a constant."):
							else
								tmppset:=[op(tmppset), DP]:
								userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"cell", cellIndex, "is zero-dim a delineating polynomial was used."):
								userinfo(4, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"the delineating polynomial used is", DP):
							fi:
						else 
							userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "there is nullification on cell", cellIndex, "which has dim>0. Hence the CAD cannot be guaranteed sign-invariant."):
							userinfo(4, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "the polynomial", op(k,pset[i]), "is nullified by", alpha, "which is a cell of dim>0."):
							if failure='err' then
								error("The input is not well-oriented (there is nullification on cell %1).  The output cannot be guaranteed correct.", cellIndex):
							elif failure='giveFAIL' then
								return(FAIL):
							else
								WARNING("The input is not well-oriented (there is nullification on cell %1).  The output cannot be guaranteed correct.", cellIndex):
								tmppset:=[op(tmppset)]:
							fi:
						fi:
					else 
						tmppset:=[op(tmppset)]: 
						userinfo(3, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"the nullified polynomial on cell", cellIndex, "was discarded since we are using Collins lifting."):
					fi:
				fi: 
			else tmppset:=[op(tmppset),op(k,pset[i])]:
			fi:
		od: 
		tmpcad[j]:=op(PCAD_GenerateStack( cell, tmppset, R[i], out )):
	od: 
	cad[i]:=[seq( tmpcad[k], k=1..nops(cad[i-1]) )]: 
	cad[i]:=map(X->op(1,X),cad[i]): 
	if cad[i]<>sort(cad[i], 'PCAD_SortCadCell') then error("the cells are not in order. This is an unexpected error message - please report to M.England@bath.ac.uk"); fi: 
	userinfo(2, {'ProjectionCAD', 'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"produced CAD of", R[i]['variables'], "-space with", nops(cad[i]), "cells"):
	if retcad=i then 
		return(cad[i]):
	fi:
od:
if OMswitch=0 then 
	kernelopts('opaquemodules'=true): 
fi:
if out='list' or out='listwithrep' then
	return( cad[dim] ):
elif out='rootof' then 
	tmpcad:=map(X -> X[2],cad[dim]):
	tmpcad:=map(Y-> map(X -> `if`(op(0,X)=And, op(X), X), Y), tmpcad);
	return(tmpcad):
elif out='piecewise' then
	return( PCAD_LWRCADtoPWCAD( cad[dim] ) ): 
fi:
end proc:

#####################################################################

PCAD_LWRCADtoPWCAD:=proc( LCAD::list ) :: 'list':
local dim, Tree, Leaf, i, NumRoots, SubCAD, Rep:
description "PCAD_LWRCADtoPWCAD: Convert a Cad in listwithrep output format to a piecewise CAD.":

dim:=convert(map(X->nops(X[1]),LCAD),set):
if nops(dim)<>1 then
	error("cells have indices of different length. This is an unexpected error message - please report to M.England@bath.ac.uk"):
else
	dim:=op(dim):
fi:
if dim=0 then
	if nops(LCAD)<>1 or LCAD[1][1]<>[] or LCAD[1][2]<>[] then
		error("converting to piecewise did not go as expected. This is an unexpected error message - please report to M.England@bath.ac.uk"):
	fi: 
	Leaf:=LCAD[1][3]:
	return(Leaf):
fi:
NumRoots:=max(map(X->X[1,1],LCAD)):
Tree:=[]:
for i from 1 to NumRoots do 
	SubCAD:=select(X->X[1][1]=i, LCAD):
	Rep:=convert(map(X->X[2,1],SubCAD),set):
	if nops(Rep)<>1 then
		error("cells on same branch have different representations. This is an unexpected error message - please report to M.England@bath.ac.uk"):
	else
		Rep:=op(Rep):
	fi:
	SubCAD:=map(X->[X[1][2..dim], X[2][2..dim], X[3]], SubCAD):
	Tree:=[op(Tree), Rep, PCAD_LWRCADtoPWCAD( SubCAD )]: 
od:
piecewise(op(Tree)):
end proc:

#####################################################################
#  Section 6 - TTICAD and ECCAD
#####################################################################

TTI_TTIProjectionOperator := proc( E::table, A::table, t::posint, mvar::symbol, lvars::list, $ ) :: set:
local i, j, BigE, RESX, P, BigP:
description "TTI_TTIProjectionOperator: Computes the application of the TTICAD Projection operator.",
			"Input: Tables E and A, an integer t, the main variable for the application and a list of any other variables.  The integr t is the number of QFFs in Phi.",  
			"The tables have entries for 1..t with each entry a set of polynomials.",
			"Output: The set of polynomials obtained by applying the TTI projection operator.":

BigE := map(X -> op(X), {entries(E, 'nolist')}):
RESX := {}:
for i from 1 to nops(BigE) do 
	for j from i+1 to nops(BigE) do:
		RESX:={op(RESX), resultant(op(i,BigE), op(j,BigE), mvar)}:
	od:
od:
P := table():
for i from 1 to t do 
	P[i]:=TTI_ECReducedProjectionOperator( E[i], A[i], mvar, lvars):
od:
BigP := `union`(seq(P[i], i=1..t)):
BigP := BigP union RESX:
BigP := remove(X -> X::'constant', BigP):
RETURN(BigP):
end proc:

#####################################################################

TTI_ECReducedProjectionOperator := proc( E::set, A::set, mvar::symbol, lvars::list, $ ) :: set:
local j, k, numE, numA, P:
description "TTI_ECReducedProjectionOperator: Computes the application of the reduced projection operator for equational constraints.",
			"Input: Sets of polynomials E and A, the main variable for the application and a list of any other variables. ",  
			"Output: The set of polynomials obtained by applying the reduced projection operator.":

numE:=nops(E):
numA:=nops(A):
P := PCAD_McCallumProj(E, mvar, lvars):
for j from 1 to numE do 
	for k from 1 to numA do:
		P:={op(P), resultant(op(j,E), op(k,A), mvar)}:
	od:
od:
P := remove(X->X::'constant', P):
end proc:

#####################################################################

TTI_ECCAD:=proc( IN::list, vars::list(symbol), failure, LiftAll::boolean, output::symbol, retcad::nonnegint, $ )
local f, gs, i, n, mvar, lvars, E, F, A, B, C, P, lowercad, cad, cell, FinalCAD:
description "TTI_ECCAD: Compute the CAD using equational constraint at first level.",
			"Input: f, a polynomial representing an equational constraint and gs, a list of polynomials representing the other constraints.",
			"Also a list of ordered variables and optionally, an ouput choice and retcad level."
			"Output: A CAD which is sign invariant for f and also for g on those cells where f=0.":

f:=IN[1]:
gs:=IN[2]:
n := nops(vars):
mvar := vars[1]:
lvars := remove(X->X=mvar,vars):
E := {f}:
F := PCAD_PrimSet(E, mvar):
F := PCAD_SFBasis(F, mvar):
if n=1 then 
	return( TTI_UnivariateCase(F, mvar, output) ):
fi:
A := {f, op(gs)}:
B := PCAD_PrimSet(A, mvar):
C := PCAD_ContSet(A, mvar):
P := TTI_ECReducedProjectionOperator( F, B, mvar, lvars):
userinfo(4, {'ProjectionCAD'},"applying the reduced projection operator gave the following projection polynomials:", P):
lowercad := TTI_LowerDimCAD(P union  C, lvars, output, retcad, failure): 
if lowercad=FAIL then 
	print("The reduced projection set is not well oriented.  Hence CADW fails and the ECCAD cannot be produced"):
	return(FAIL):
fi:
if retcad=n or retcad=0 then
	cad:=table():
	for i from 1 to nops(lowercad) do
		cell := lowercad[i]:
		cad[i] := TTI_TTIGenerateStack( cell, 1, mvar, lvars, table([(1)=gs]), table([(1)=B]), table([(1)=E]), LiftAll, output, failure):
		if cad[i]=FAIL then
			print("The equational constraint is nullified on a cell of positive dimension and the excluded polynomials are not constant.  Hence this implementation of ECCAD fails."):
			return(FAIL):
		fi:
	od:
	FinalCAD := [entries(cad, 'nolist')]:
	FinalCAD := map(X->op(X), FinalCAD):
	userinfo(2, {'ProjectionCAD'},"produced CAD of", vars, "-space with", nops(FinalCAD), "cells"):
else 
	FinalCAD := lowercad:
fi:
if output='list' or output='listwithrep' then 
	return(FinalCAD):
elif output='rootof' then 
	return( map(X->X[2], FinalCAD) ):
elif output='piecewise' then 
	return( PCAD_LWRCADtoPWCAD(FinalCAD) ):
else 
	error("unexpected output format.  This is an unexpected error, please report to M.England@bath.ac.uk"):
fi:
end proc:

#####################################################################

TTI_ECPP:=proc( IN::list, vars::list(symbol), $ ) :: set:
local f, gs, inpols, n, mvar, lvars, E, F, A, B, C, topP, lowP, PP:
description "TTI_ECPP: Compute the set of polynomials used by ECCAD for pojection.",
			"Input: f, a polynomial representing an equational constraint and gs, a list of polynomials representing the other constraints.",
			"Also a list of ordered variables.",
			"Output: The set of projection polynomials.":

f := IN[1]:
gs := IN[2]:
n := nops(vars):
inpols := {f,op(gs)}:
mvar := vars[1]:
lvars := remove(X->X=mvar,vars):
E := {f}:
F := PCAD_PrimSet(E, mvar):
F := PCAD_SFBasis(F, mvar):
if n=1 then
	PP:=E:
else
	A := {f, op(gs)}:
	B := PCAD_PrimSet(A, mvar):
	C := PCAD_ContSet(A, mvar):
	topP := TTI_ECReducedProjectionOperator( F, B, mvar, lvars):
	lowP := PCAD_McCProjPolys(topP union C, lvars):
	PP := topP union lowP union inpols:
fi:
PCAD_SetFactors(PP):
end proc:

#####################################################################

TTI_TTICAD:=proc( PHI::list, vars::list, failure, LiftAll::boolean, output::symbol, retcad::nonnegint, $ )
local i, t, n, mvar, lvars, phi, f, gs, E, F, A, B, C, BigF, BigC, P, lowercad, cad, cell, FinalCAD:
description "TTI_TTICAD: Compute the TTICAD for a list of QFFs.",
			"Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial (equational constraint) and the second a list of polynomials (from the other constraints).",
			"Also a list of ordered variables and optionally, an ouput choice and retcad level."
			"Output: A truth-table invariant CAD for PHI or failure if the polynomials are not well oriented.":

t := nops(PHI):
n := nops(vars):
mvar := vars[1]:
lvars := remove(X->X=mvar,vars):
phi := table():
f := table():
gs := table():
E := table():
F := table():
for i from 1 to t do 
	phi[i] := PHI[i]:
	f[i] := expand(phi[i][1]):
	gs[i] := expand(phi[i][2]):
	if f[i]=false then
		E[i] := {op(gs[i])}:
	else
		E[i] := {f[i]}:
	fi:
	F[i] := PCAD_PrimSet(E[i], mvar):
	F[i] := PCAD_SFBasis(F[i], mvar):
od:
BigF := `union`(seq(F[i], i=1..t)):
if n=1 then 
	return( TTI_UnivariateCase(BigF, mvar, output) ):
fi:
A := table():
B := table():
C := table():
for i from 1 to t do 
	A[i] := {f[i], op(gs[i])}:
	B[i] := PCAD_PrimSet(A[i], mvar):
	C[i] := PCAD_ContSet(A[i], mvar):
od:
BigC := `union`(seq(C[i], i=1..t)):
P := TTI_TTIProjectionOperator(F, B, t, mvar, lvars):
P := P union  BigC:
userinfo(4, {'ProjectionCAD'},"applying the reduced projection operator gave the following projection polynomials:", P):
lowercad := TTI_LowerDimCAD(P, lvars, output, retcad, failure): 
if lowercad=FAIL then 
	WARNING("The reduced projection set is not well oriented.  Hence CADW fails and the TTICAD cannot be produced"):
	return(FAIL):
fi:
if retcad=0 or retcad=n then
	cad:=table():  
	for i from 1 to nops(lowercad) do
		cell := lowercad[i]:
		cad[i] := TTI_TTIGenerateStack( cell, t, mvar, lvars, gs, B, E, LiftAll, output, failure):
		if cad[i]=FAIL then
			WARNING("An equational constraint is nullified on a cell of positive dimension and the excluded polynomials are not constant.  Hence this implementation of TTICAD fails."):
			return(FAIL):
		fi:
	od:
	FinalCAD := [entries(cad, 'nolist')]:
	FinalCAD := map(X->op(X), FinalCAD):
	userinfo(2, {'ProjectionCAD'},"produced CAD of", vars, "-space with", nops(FinalCAD), "cells"):
else
	FinalCAD := lowercad:
fi:
if output='list' or output='listwithrep' then 
	return(FinalCAD):
elif output='rootof' then 
	return( map(X->X[2], FinalCAD) ):
elif output='piecewise' then 
	return( PCAD_LWRCADtoPWCAD(FinalCAD) ):
else 
	error("unexpected output format.  This is an unexpected error, please report to M.England@bath.ac.uk"):
fi:
end proc:

#####################################################################

TTI_TTIPP:=proc( PHI::list, vars::list, $ ) :: set:
local i, t, n, mvar, lvars, phi, f, gs, E, F, A, B, C, BigF, BigC, BigA, topP, lowP, PP:
description "TTI_TTIPP: Compute the projection polynomials used in TTICAD for a list of QFFs.",
			"Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial (equational constraint) (or false if there is none) and the second a list of polynomials (from the other constraints).",
			"Also a list of ordered variables."
			"Output: The set of projection polynomials.":

t := nops(PHI):
n := nops(vars):
mvar := vars[1]:
lvars := remove(X->X=mvar,vars):
phi := table():
f := table():
gs := table():
E := table():
F := table():
for i from 1 to t do 
	phi[i] := PHI[i]:
	f[i] := expand(phi[i][1]):
	gs[i] := expand(phi[i][2]):
	if f[i]=false then
		E[i] := {op(gs[i])}:
	else
		E[i] := {f[i]}:
	fi:
	F[i] := PCAD_PrimSet(E[i], mvar):
	F[i] := PCAD_SFBasis(F[i], mvar):
od:
BigF := `union`(seq(F[i], i=1..t)):
if n=1 then 
	PP:=BigF:
else
	A := table():
	B := table():
	C := table():
	for i from 1 to t do 
		A[i] := {f[i], op(gs[i])}:
		B[i] := PCAD_PrimSet(A[i], mvar):
		C[i] := PCAD_ContSet(A[i], mvar):
	od:
	BigA := `union`(seq(A[i], i=1..t)):
	BigC := `union`(seq(C[i], i=1..t)):
	topP := TTI_TTIProjectionOperator(F, B, t, mvar, lvars):
	lowP := PCAD_McCProjPolys(topP union BigC, lvars):
	PP:= topP union lowP union BigA
fi:
PCAD_SetFactors(PP):
end proc:

#####################################################################

TTI_UnivariateCase := proc( F::set(polynom), var::symbol, output::symbol, $ )
local cad:
description "TTI_UnivariateCase: Dealing with the case where there is only one variable.":

if output='list' then 
	cad := RegularChains:-TRDgenerate_stack( [], convert(F,'list'), RegularChains:-PolynomialRing([var]) ):
	cad:=map(X->X[1], cad):
	return(cad):
else
	cad := RegularChains:-TRDgenerate_stack_general( [], convert(F,'list'), RegularChains:-PolynomialRing([var]) ):
fi:
cad:=map(X->X[1], cad):
if output='listwithrep' then 
	return(cad):
elif output='piecewise' then 
	return( PCAD_LWRCADtoPWCAD(cad) ):
elif output='rootof' then 
	return( map(X->X[2], cad) ):
else error("unexpected output format.  This is an unexpected error, please report to M.England@bath.ac.uk"):
fi:
end proc:

#####################################################################

TTI_LowerDimCAD := proc( P::set(polynom), lvars::list(symbol), output::symbol, retcad::nonnegint, failure, $ )
local pset, cad:
description "TTI_LowerDimCAD: Constructing the lower dimension CAD using CADW via PCAD.":

pset := PCAD_McCProjPolys(P, lvars):
userinfo(2, {'ProjectionCAD'},nops(pset), "projection factors were calculated for the lower dimensional CAD."):
userinfo(4, {'ProjectionCAD'},"the factors are:", pset):
if output='list' then 
  cad:=PCAD_ProjCADLift( pset, lvars, 'McCallum', 'OI', retcad, 'list', failure ):
else
  cad:=PCAD_ProjCADLift( pset, lvars, 'McCallum', 'OI', retcad, 'listwithrep', failure ):
fi:
cad:
end proc:

#####################################################################

TTI_TTIGenerateStack := proc( cell::list, t::posint, mvar::symbol, lvars::list(symbol), gs::table, B::table, E::table, LiftAll::boolean, output::symbol, failure, $ ) :: list:
local j, vars, R, Stack, SPrc, SPcube, alpha, beta, falpha, L, cellIndex, excl, entry, nullification, extra:
description "TTI_TTIGenerateStack: Creating the lifting set and then stack over a cell in the lower dimensional CAD.":

cellIndex := cell[1]:
vars:=[mvar, op(lvars)]:
R := RegularChains:-PolynomialRing(vars):
SPrc := cell[-1][1]:
SPcube := cell[-1][2]:
alpha := PCAD_SPtoRootOf( SPrc, SPcube, RegularChains:-PolynomialRing(lvars) ):
beta := []:
for j from 1 to nops(cellIndex) do 
	if cellIndex[j]::'even' then 
		beta:=[op(beta), alpha[j]]:
	fi:
od:
L:={}:
for j from 1 to t do
	nullification := false:
	for entry in E[j] do 
		falpha := subs(alpha, entry):
		if is(falpha=0)=true then 
			userinfo(3, {'ProjectionCAD', 'TTICAD'},"the equational constraint factor", entry, " is nullified on the cell ", cell):
			nullification := true:
		fi:
	od:
	if nullification=true then 
		if PCAD_IsCellZeroDim(cellIndex)=true then
			L := L union B[j]:
			L := remove(X -> evalb(is(subs(alpha,X)=0)), L):
			userinfo(3, {'ProjectionCAD', 'TTICAD'},"the cell is zero-dimensional so we can continue by expanding the lifting set on this cell."):
		else
			excl := TTI_ExclProj_phi(gs[j], vars):
			excl := remove(X->X::'constant', subs(beta, excl)):
			userinfo(4, {'ProjectionCAD', 'TTICAD'},"the excluded polynomials are ", excl, " on the cell."):
			if excl = {} then 
				L := L union B[j]:
				L := remove(X -> evalb(is(subs(alpha,X)=0)), L):
				userinfo(3, {'ProjectionCAD', 'TTICAD'},"the relevant excluded polynomials are constants, hence the non-equational constraints are delineable on the cell and we can continue by expanding the lifting set on this cell."):
			else
				userinfo(3, {'ProjectionCAD', 'TTICAD'},"the cell has positive dimension and the relevant excluded polynomials are non-constants.  Hence the algorithm fails."):
				if failure='err' then
					error("The input is not well-oriented (there is nullification on cell %1). The output cannot be guaranteed correct.", cellIndex):
				elif failure='giveFAIL' then
					return(FAIL):
				else
					WARNING("The input is not well-oriented (there is nullification on cell %1). The output cannot be guaranteed correct.", cellIndex):
					userinfo(3, {'ProjectionCAD', 'TTICAD'},"we continue, expanding the lifting set on this cell."):
					L := L union B[j]:
					L := remove(X -> evalb(is(subs(alpha,X)=0)), L):
				fi:
			fi:
		fi:
	else
		if LiftAll=true then 
			userinfo(3, {'ProjectionCAD', 'TTICAD'},"we only need to lift with respect to the ECs, but since the LiftAll option was given we expand the lifting set."):
			# In this case we include all the B[j] in the lifting set.  This simulates QEPCAD.  It will result in extra cells making the gs sign-invariant.  
				#McCallum's theory allows us to conclude the gs are sign-invariant in SECTIONS OF F, not elsewehere.  In particular they are not delineable elsewhere. Thus the output is only sign-invariant for the sample point, not neccessarily the cell.
			extra := convert(B[j], 'list'):
			# Also, in this case we need to ensure the extra ones are coprime at squarefree over the sample point (it is not enough to do it for the real point and use delineability to argue it works for the rest of the cell).
				# We do that extra work now.
			extra := remove(X -> eval(subs(alpha,X))::'constant', extra):
			extra := PCAD_MakeCoprimeOverPoint( SPrc, SPrc, SPcube, extra, R ): 
			extra := PCAD_MakeSquareFreeOverPoint( SPrc, [], SPcube, extra, R ):
			# We then let these extra polynomials join L.
			L := L union convert(extra, 'set'):
		else
			L := L union E[j]:
		fi:
	fi:
od:
L := PCAD_SetFactors(L):
L := convert(L, 'list'):
if output=list then
	Stack := PCAD_GenerateStack(cell, L, R, 'list'):
else 
	Stack := PCAD_GenerateStack(cell, L, R, 'listwithrep'):
fi:
Stack := map(X->X[1],Stack):
end proc:

#####################################################################

TTI_ResCADSet := proc( PHI::list, vars::list(symbol), $ ) :: list:
local i, j, k, t, n, mvar, phi, f, gs, E, F, A, B, C, R, numE, numA, BigE: 
description "TTI_ResCADSet: Compute the ResCAD set.",
			"Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial (equational constraint) and the second a list of polynomials (from the other constraints).",
			"Also a list of ordered variables."
			"Output: The ResCAD set.":

t := nops(PHI):
n := nops(vars):
if n=1 then
	WARNING("There is only one variable and hence no projection and no ResCADSet."):
	return([]):
fi:
mvar := vars[1]:
phi := table():
f := table():
gs := table():
A := table():
B := table():
C := table():
E := table():
F := table():
for i from 1 to t do 
	phi[i] := PHI[i]:
	f[i] := expand(phi[i][1]):
	gs[i] := expand(phi[i][2]):
	if f[i]=false then
		E[i] := {op(gs[i])}:
		A[i] := {op(gs[i])}:
	else
		E[i] := {f[i]}:
		A[i] := {f[i], op(gs[i])}:
	fi:
	if has(map(X->has(X, mvar), E[i]), false) then
		WARNING("An equational constraint will certainly be nullified since %1 does not contain the main variable %2.", f[i], mvar):
	fi:
	F[i] := PCAD_PrimSet(E[i], mvar):
od:
BigE := `union`(seq(E[i], i=1..t)):
for i from 1 to t do 
	B[i] := PCAD_PrimSet(A[i], mvar):
	C[i] := PCAD_ContSet(A[i], mvar):
od:
R:={}:
for i from 1 to t do 
	numE := nops(E[i]):
	numA := nops(A[i]):
	for j from 1 to numE do 
		for k from 1 to numA do:
			R := {op(R), resultant(op(j, E[i]), op(k, A[i]), mvar)}:
		od:
	od:
od:
R := R union BigE:
PCAD_RemoveConstantMultiples( remove(X->X::'constant',R) ):
end proc:

#####################################################################

TTI_ResCAD := proc( PHI::list, vars::list(symbol), output, $ ) :: list:
local i, t, n, mvar, phi, f, gs, E, F, R, BigF, pset, cad: 
description "TTI_ResCAD: Compute a TTICAD for a list of QFFs using the RESCAD method.",
			"Input: PHI, a list of lists phi[i], each of which represent a QFF.  Each phi[i] has two arguments; the first a polynomial (equational constraint) and the second a list of polynomials (from the other constraints).",
			"Also a list of ordered variables and optionally, an ouput choice."
			"Output: Failure if the reduced projection set is not well oriented.  Otherwise a CAD for PHI which is guarenteed truth table invariant provided that no equational constraint is nullified.":

n := nops(vars):
if n=1 then 
	t := nops(PHI):
	mvar := vars[1]:
	phi := table():
	f := table():
	gs := table():
	E := table():
	F := table():
	for i from 1 to t do 
		phi[i] := PHI[i]:
		f[i] := expand(phi[i][1]):
		if has(f[i], mvar)=false then
			WARNING("An equational constraint will certainly be nullified since %1 does not contain the main variable %2.", f[i], mvar):
		fi:
		gs[i] := expand(phi[i][2]):
		E[i] := {f[i]}:
		F[i] := PCAD_PrimSet(E[i], mvar):
	od:
	BigF := `union`(seq(F[i], i=1..t)):
	return( TTI_UnivariateCase(BigF, mvar, output) ):
fi:
R := TTI_ResCADSet(PHI, vars):
pset := PCAD_McCProjPolys(R, vars):
cad := PCAD_ProjCADLift(pset, vars, 'McCallum', 'SI', 0, output, 'giveFAIL' ):
if cad=FAIL then
	print("The reduced projection set is not well oriented.  Hence CADW fails and the ResCAD cannot be produced"):
	return(FAIL):
fi:
cad:
end proc:

#####################################################################

TTI_ExclProj_phi := proc(gs::list, vars::list(symbol)) :: set:
local mvar, excl, ent:
description "TTI_ExclProj_phi: Compute the set of polynomials excluded from the full projection set for a QFF by the redcued operator.",
			"Input: gs the set of non equational constraints for a QFF and the list of variables.",
			"Output: Set of elculded polynomials.":

mvar := vars[1]:
excl := {}:
for ent in gs do
	excl := {op(excl), discrim(ent, mvar), coeffs(expand(ent), mvar) }:
od:
excl := remove(X->X::'constant', excl):
excl := PCAD_RemoveConstantMultiples(excl):
PCAD_SFBasis(excl, vars[2]);
end proc:

#####################################################################
#  Section 7 - Formulation Commands
#####################################################################

Formulation_VarOrdFree :=  proc( vars::list(symbol), $ ) :: list:
description "Formulation_VarOrdFree: Finds all the possible variable orderings for a list of variables.",
			"Input: A list of variables.",
			"Ouput: A list of all posible acceptable variable orderings (each a list of variables).":

eval('combinat'['permute'](vars));
end proc:

#####################################################################

Formulation_VarOrdRestricted :=  proc( vars::list(list(symbol)), $ ) :: list:
local i, nBlocks, Parts, Poss, nPoss, num:
description "Formulation_VarOrdFree: Finds all the possible variable orderings for a list of variables in blocks.",
			"Input: A list of lists of variables.",
			"Ouput: A list of all posible acceptable variable orderings (each a list of variables).":

nBlocks:=nops(vars):
Parts:=table():
for i from 1 to nBlocks do
	Parts[i] := Formulation_VarOrdFree( vars[i] ):
od:
Poss:=Parts[1]:
num:=2:
while num<nBlocks+1 do
	nPoss:=table(): 
	for i from 1 to nops(Parts[num]) do
		 nPoss[i] := map(X -> [op(X), op(Parts[num][i])], Poss):
	od:
	Poss := map(X->op(X), [entries(nPoss, 'nolist')]):
	num:=num+1:
od:
Poss:
end proc:

#####################################################################

Formulations_RGS := proc(N::posint, i::posint) :: list:
local j, s, f, L, next_str:
description "Formulations_RGS: Creates the restricted growth strings for given input.",
			"Input: N, the length of the strings and i, the maximum permitted block size.",
			"Ouput: List of possible restricted growth strings.  When i=1 we get the rgs' for possible partitions of a set of size N.":

next_str := proc() 
local k, j, sk, m1, mp:
# Returns index of first changed element in s[] and zero if current string is the last
k := N:
do
	k := k-1:
	if k = 0 then 
		return(0):
	fi:
	sk := s[k] + 1;
	m1 := f[k - 1];
	mp := m1 + i;
	if sk > mp then
		s[k] := 0;
	else
		s[k] := sk:
		if sk = mp then 
			m1 := m1 + i: 
		fi:
		for j from k to N-1 do 
			f[j] := m1:
		od:
		return(k):
	fi
od
end proc:
s := array(0..N-1);
f := array(0..N-1);
for j from 0 to N-1 do 
	s[j] := 0:
od:
for j from 0 to N-1 do 
	f[j] := 0:
od:
L := NULL:
do
	L := L, convert(s,list);
	if next_str() = 0 then 
		break 
	fi:
od:
[L]:
end proc:

#####################################################################

Formulations_Classification := proc(s::list) :: list: 
local L, i, m:
description "Formulations_Classification: Turns a restricted growth string (in the case i=1) to a set partition.",
			"Input: s, a restricted growth string.",
			"Ouput: The set partition corresponding to s.":

m := 1 + max(op(s));
L:=table(): 
for i from 1 to m do 
	L[i]:=[]: 
od: 
L:=[entries(L, 'nolist')]:
for i from 1 to nops(s) do
	L[s[i]+1] := [op(L[s[i]+1]),i];
od;
L:
end proc:

#####################################################################

Formulations_SetPartitions := proc(n::posint) :: list: 
local P, S, s, i:
description "Formulations_SetPartitions: Create the possible set partitions of an integer.",
			"Input: n, a positive integer.",
			"Ouput: The list of all set partitions of n.":
P := table(): 
S := Formulations_RGS(n,1);
for i from 1 to nops(S) do 
  s:=S[i]:
  P[i]:=Formulations_Classification(s):
od:
[entries(P, 'nolist')]:
end proc:

#####################################################################

Formulations_FunctionPartitions := proc(L::list) :: list:
local SP, entry, sets, i, count, s:
description "Formulations_FunctionPartitions: Create the possible set partitions of a list of functions.",
			"Input: L, a list of functions.",
			"Ouput: The list of all set partitions (stored as lists) of the functions.":

SP := Formulations_SetPartitions(nops(L)):
sets := table(): count:=1:
for entry in SP do
	for i from 1 to nops(entry) do
		s[i]:=L[op(i,entry)]:
	od:
	sets[count] := [seq(s[i], i=1..nops(entry))]:
	count := count+1:
od:
[entries(sets, 'nolist')]:
end proc:

#####################################################################

Formulations_ECsForTTI := proc(L::list) :: list:
local fparts, ECparts, entry, i, SS, entrysets, j, newentrysets, k, l, ent, count, COUNT;
description "Formulations_ECsForTTI: Create the list of partitions of equational constraints for use in TTICAD.",
			"Input: L, a list of polynomials (assumed equational constraints).",
			"Ouput: The list of all possible partitions (stored as lists) in which the first entry is a designated equational constraint.":

fparts := Formulations_FunctionPartitions(L):
ECparts := table(): 
COUNT:=1:
for entry in fparts do
	for i from 1 to nops(entry) do
		SS[i] := [seq( [entry[i][j], op(remove(has,entry[i],entry[i][j]))], j=1..nops(entry[i]) )]:
	od;
	entrysets := [[]]:
	for j from 1 to nops(entry) do
		newentrysets := table(): 
		count:=1:
		for k from 1 to nops(SS[j]) do
			for l from 1 to nops(entrysets) do
				ent := [op(entrysets[l]), SS[j][k]]:
				newentrysets[count] := ent:
				count := count+1:
			od:
		od:
	entrysets := [entries(newentrysets, 'nolist')]:
	od:
	ECparts[COUNT] := [op(entrysets)]:
	COUNT := COUNT+1:
od:
map(X->op(X), [entries(ECparts, 'nolist')]):
end proc:

#####################################################################

Formulations_ECCAD := proc( IN::list ) :: list:
local i, ecs, necs:
description "Formulation_ECCAD: Create the list of choices for entry into ECCAD.",
			"Input: Two lists of polynomials, the first assumed equational constraints.",
			"Ouput: A list of possible formulations for ECCAD.  Each list has a polynomial (designated EC) and a list of other constraints.":

ecs := IN[1]:
necs := IN[2]:
[seq( [ecs[i],[op(remove(X->X=ecs[i],ecs)),op(necs)]], i=1..nops(ecs))]:
end proc:

#####################################################################

Formulations_QFF := proc( IN::list(list(polynom)), $ ) :: list:
local F, G, count, Parts, i, j, k, tempParts:
description "Formulations_Claise: Create the list of possible formulations for a TTICAD QFF.",
			"Input: Two lists of polynomials, the first assumed the equational constraints.",
			"Ouput: The list of all possible formulations for TTICAD (inc splitting).":

F := IN[1]:
G := IN[2]:
Parts := Formulations_ECsForTTI(F):
for i from 1 to nops(G) do
	tempParts := table(): 
	count := 1:
	for j from 1 to nops(Parts) do
		for k from 1 to nops(Parts[j]) do
			tempParts[count] := [op(1..k-1, Parts[j]), [op(Parts[j][k]), G[i]], op(k+1..nops(Parts[j]), Parts[j])]:
			count := count+1:
		end do:
	end do:
	Parts := [entries(tempParts, 'nolist')]:
od:
map(X-> map(Y-> [Y[1],Y[2..nops(Y)]] ,X), Parts):
end proc:

#####################################################################

Formulations_TTICAD := proc( LL::list(list(list(polynom))), splitting::boolean, $ ) :: list:
local i, N, F, G, Parts, Poss, nPoss, num:
description "Formulations_TTICAD: Create the list of possible formulations for TTICAD.",
			"Input: A list of lists each representing a QFF.  Each QFF has two lists of polynomials, the first assumed the equational constraints.  Optionally, the choice of whether to consider splitting.",
					"Note: it is assumed that the QFFs come from a disjunctive normal form.  I.e. trivial merging is not possible.",
			"Ouput: The list of all possible formulations for TTICAD, by default including splitting.":

N := nops(LL):
if splitting=false then
	for i from 1 to N do
		F:=LL[i][1]: 
		G:=LL[i][2]:
		Parts[i] := map(X->[X], Formulations_ECCAD([F,G])):
	od:
else 
	for i from 1 to N do
		F:=LL[i][1]: 
		G:=LL[i][2]:
		Parts[i] := Formulations_QFF( [expand(F), expand(G)] ):
	od:
fi:
Poss:=Parts[1]:
num:=2:
while num<N+1 do
	nPoss:=table(): 
	for i from 1 to nops(Parts[num]) do
		 nPoss[i]:=map(X->[op(X),op(Parts[num][i])], Poss):
	od:
	Poss:=map(X->op(X), [entries(nPoss, 'nolist')]):
	num:=num+1:
od:
Poss:
end proc:

#####################################################################
#  Section 8 - Heuristics
#####################################################################

Heuristics_sotdP := proc(g::polynom) :: posint:
local G:
description "Heuristics_sotdP: Calculate the sotd (sum of total degree) of a polynomial.":

G:=expand(g):
if G::`+` then 
	return( map(degree, expand(g)) );
else 
	return( degree(g) ):
fi:
end proc:

#####################################################################

Heuristics_sotdL := proc(IN) :: posint:
local L, f:
description "Heuristics_sotdL: Calculate the sotd (sum of total degree) of a list / set of polynomials.":

if IN::'set' then 
	L := convert(IN, 'list'):
elif IN::'list' then
	L := IN:
else
	error("invalid input:  sotdL expects its input to be a list or set of polynomials.  This is an unexpected erroor - please report to M.England@bath.ac.uk"):
fi:
add( Heuristics_sotdP(f), f in L):
end proc:

#####################################################################

Heuristics_ndrr := proc(INN, var::symbol, method::symbol, naivedeg::posint, Quiet::boolean, TimeO::nonnegint, $) :: posint:
local ans, i, othervars, Uni, pol, deg, SFUni:
description "Heuristics_ndrr: Calculate the ndrr (number of distinct real roots) of a polynomial or list / set of polynomials with respect to a variable.  Note: if multivariate polynomials are included then they are ignored.",
			"Input: Either a polynomial of list / set of polynomials, and a variable.",
			"Ouput: The ndrr.":

if TimeO<>0 then
	try
		ans := timelimit(TimeO, Heuristics_ndrr(INN, var, method, naivedeg, Quiet, 0)):
	catch "time expired":
		ans := infinity:
	end try:
	return(ans):
fi:
if INN::'polynom' then 
	if indets(INN, 'name') <> {var} then
		return(0):
	else 
		return( sturm(INN, var, -infinity, infinity) ):
	fi:
elif INN::'set' or INN::'list' then 
	othervars := remove(X->X=var,indets(INN,'name')):
	Uni := convert(remove(has,INN,othervars),'set'):
	if Uni={} then
		return(0):
	fi:
	if method='naive' then
		return( add( sturm(Uni[i], var, -infinity, infinity), i=1..nops(Uni) ) ):
	elif method='SFnaive' then 
		SFUni := PCAD_SFBasis(Uni, var):
		return( add( sturm(SFUni[i], var, -infinity, infinity), i=1..nops(Uni) ) ):
	else 
		pol := mul(Uni[i], i=1..nops(Uni)):
		deg := degree(pol):
		if method='full' then
			return( sturm( pol, var, -infinity, infinity) ):
		elif method='vary' then
			if deg < naivedeg then 
				return( sturm( pol, var, -infinity, infinity) ):
			else 
				if Quiet=false then
					WARNING("The method is %1 and degree of the product of the polynomials is greater that naivedeg=%2.  Hence the ndrr for the individual polynomials will been returned, which is an overestimate.  To avoid this use the optional argument method=full or naivedeg=int where int>80", method, naivedeg): 
				fi:
				return( add( sturm(Uni[i], var, -infinity, infinity), i=1..nops(Uni) ) ):
			fi:
		else 
			error("this is an unexpected error.  Please report to M.England@bath.ac.uk"):
		fi:
	fi:
else
	error("this is an unexpected error.  Please report to M.England@bath.ac.uk"):
fi:
end proc:

#####################################################################

Heuristics_BHFeatures := proc(INN, var::symbol, $) :: list:
local i, polL, degs, varonly, nterms, BH1, BH2, BH3, f:
description "Heuristics_BHFeatures: Calculate the three features used for Brown's Heuristic.",
			"Input: A list of polynomials, and a variable.",
			"Ouput: A list of three positive integers: (1) The highest degree of the variable occuring, (2) The highest total degree of a term the variable appears in, (3) The number of terms with that variable.":

if INN::'list' then
	polL := expand(INN):
elif INN::'set' then
	polL := convert(expand(INN), 'list'):
else
	error("First argument should be list or set of polynomials"):
fi:
degs := map(X -> degree(X, var), polL ):
BH1 := max(degs):
varonly := map(X -> select(has, X, var), polL);
degs := map(X -> degree(X), varonly ):
BH2 := max(degs):
f := proc(p): if type(p,`+`) then return(nops(p)) else return(1): fi: end proc:
nterms := map(X -> f(X), varonly):
BH3 := add( nterms[i], i=1..nops(nterms) ):
return( [BH1, BH2, BH3] );
end proc:

#####################################################################

Heuristics_PickSuggestions := proc( Poss::list, Hvalues::list, {pick::symbol:='minimise'}, $) :: list:
local bestH, Suggestions, count, i:
description "Heuristics_PickSuggestions: Selects suggested formulations based on heuristic values.",
			"Input: Two lists, the first a list of possible formulations, the second a list of heuristic values.  It is assumed we are picking the formulations to minimise the heuristic, if maximise then specify as an option.",
			"Ouput: A list of those formulations with the min / max heuristic value.":

if nops(Poss)<>nops(Hvalues) then 
	error("unexpected error - please report to M.England@bath.ac.uk"): 
fi:
if 'pick'='minimise' then
	bestH := min(Hvalues):
elif 'pick'='maximise' then
	bestH := max(Hvalues):
else
	error("the optional input pick should be minimise or maximise"):
fi:
Suggestions := table():
count:=1:
for i from 1 to nops(Poss) do 
	if Hvalues[i]=bestH then
		Suggestions[count] := Poss[i]:
		count:=count+1:
	fi:
od:
map(X->rhs(X), sort([entries(Suggestions, 'pairs')]));
end proc:

#####################################################################

Heuristics_RunOnPPs := proc(Poss::list, PPs::table, vars::list(symbol), heuristic::symbol, ndrrOptions::list, NTimeO::nonnegint, SeeAll::boolean, Wratio::list(nonnegint), VarOrds::boolean, $) :: list(list):
local i, nPoss, Ns, Ss, Suggestions, SuggInd, Nsort, Ssort, Ws, num:
description "Heuristics_RunOnPPs: Given a table of projection polynomials and heuristic choice returns the suggestion.",
			"Input: A list of possible formulations, table of projection polynomials, choice of heuristic and whether to return one or all.",
			"Ouput: A single formulation or list of formulations suggested by the heuristic.":

nPoss := nops(Poss): 
if heuristic in {'N','SN','NS','W'} then
	if VarOrds=false then 
		Ns := [seq( Heuristics_ndrr(PPs[i], vars[-1], ndrrOptions[1], ndrrOptions[2], true, NTimeO), i=1..nPoss)]:
#		Ns := [seq( ndrr(PPs[i], vars[-1], 'method'=ndrrOptions[1], 'naivedeg'=ndrrOptions[2], 'timeout'=NTimeO), i=1..nPoss)]:
	else
		Ns := [seq( Heuristics_ndrr(PPs[i], Poss[i][-1], ndrrOptions[1], ndrrOptions[2], true, NTimeO), i=1..nPoss)]:
#		Ns := [seq( ndrr(PPs[i], Poss[i][-1], 'method'=ndrrOptions[1], 'naivedeg'=ndrrOptions[2]), i=1..nPoss)]:
	fi:
fi:
if heuristic in {'S','SN','NS','W'} then
	Ss := [seq( sotd(PPs[i]), i=1..nPoss)]:
fi:
if heuristic in {'N','NS'} then
	Suggestions:=Heuristics_PickSuggestions( Poss, Ns, 'pick'='minimise' ):
	num:=nops(Suggestions):
	userinfo(2, {'ProjectionCAD'}, "the fomulations have been compared for lowest ndrr value, identifying ", num):
	if nops(Suggestions)=1 then 
		return( op(Suggestions) ):
	elif heuristic='N' then
		return( Heuristics_Display(Suggestions, SeeAll, heuristic) ):
	else
		SuggInd := select( X -> has( map(Y->is(Y=Poss[X]),Suggestions) ,true)=true, [seq(1..nPoss)] ):
#		SuggInd := select( X -> has(Suggestions,Poss[X])=true, [seq(1..nPoss)] ):
		Ss := Ss[SuggInd]:
		Suggestions := Heuristics_PickSuggestions( Suggestions, Ss, 'pick'='minimise' ):
		num := nops(Suggestions):
		userinfo(2, {'ProjectionCAD'}, "the formulations with lowest ndrr value have had their sotd values compared for the lowest, identifying, ", num):
		return( Heuristics_Display(Suggestions, SeeAll, heuristic) ):
	fi:
elif heuristic in {'S','SN'} then
	Suggestions := Heuristics_PickSuggestions( Poss, Ss, 'pick'='minimise' ):
	num := nops(Suggestions):
	userinfo(2, {'ProjectionCAD'}, "the fomulations have been compared for lowest sotd value, identifying ", num):
	if nops(Suggestions)=1 then 
		return( op(Suggestions) ):
	elif heuristic='S' then
		return( Heuristics_Display(Suggestions, SeeAll, heuristic) ):
	else
		SuggInd := select( X -> has( map(Y->is(Y=Poss[X]),Suggestions) ,true)=true, [seq(1..nPoss)] );
#		SuggInd := select( X -> has(Suggestions,Poss[X])=true, [seq(1..nPoss)] ):
		Ns := Ns[SuggInd]:
		Suggestions:=Heuristics_PickSuggestions( Suggestions, Ns, 'pick'='minimise' ):
		num:=nops(Suggestions):
		userinfo(2, {'ProjectionCAD'}, "the formulations with lowest sotd value have had their ndrr values compared for the lowest, identifying, ", num):
		return( Heuristics_Display(Suggestions, SeeAll, heuristic) ):
	fi:
elif heuristic = 'W' then
	Ws := table():
	Nsort := sort(Ns):
	Ssort := sort(Ss):
	Ws := [seq( Heuristics_FindPos(Nsort,Ns[i])*Wratio[1] + Heuristics_FindPos(Ssort,Ss[i])*Wratio[2], i=1..nPoss)]:
	Suggestions := Heuristics_PickSuggestions( Poss, Ws, 'pick'='minimise' ):
	num:=nops(Suggestions):
	userinfo(2, {'ProjectionCAD'}, "the formulations have been ranked with respect to ndrr and sotd values with a weighted average taken according to ratio N:S = ", Wratio, ".  Picking the lowest value of this identifies ", num):
	return( Heuristics_Display(Suggestions, SeeAll, heuristic) ):
elif (heuristic='BrownFull' or heuristic='BrownBasic') then
	error("unexpected error with Brown's Heuristic - please report to M.England@bath.ac.uk"):
	# Implemented on its own.
else
	error("unexpected error - please report to M.England@bath.ac.uk"):
fi:
end proc:

#####################################################################

Heuristics_FindPos := proc( Hvalues::list(nonnegint), H::nonnegint ) :: posint:
local i:
description "Heuristics_FindPos: Find where a value is in a list.",
			"Input: A list of values and one value from that list.",
			"Ouput: The position of that value in the list.":

for i from 1 to nops(Hvalues) do
	if Hvalues[i]=H then
		return(i):
	fi:
od:
error("unexpected error - please report to M.England@bath.ac.uk"):
end proc:

#####################################################################

Heuristics_Display := proc( Suggestions::list, SeeAll::boolean, heuristic::symbol )
local N:
description "Heuristics_Display: Given the formulations suggested by a heuristic decides whether to display one or all.",
			"Input: A list of formulations, the choice of how many to see and the heuristic used.",
			"Ouput: A single formulation or a list of formulations.":

N:=nops(Suggestions):
if N=1 then
	return(op(Suggestions)):
else
	if SeeAll=true then 
		WARNING("There are %1 formulations which heuristic %2 cannot differentiate.", N, heuristic):
		return(Suggestions):
	elif SeeAll=false then
		WARNING("There are %1 formulations which heuristic %2 cannot differentiate.  Only the first has been returned.  To display them all run with option SeeAll=true.", N, heuristic):
		return( op(1,Suggestions) ):
	else 
		return( op(1,Suggestions) ):
	fi:
fi:
end proc:

#####################################################################

Heuristics_ECCAD := proc(ecs::list(polynom), necs::list(polynom), vars::list(symbol), heuristic::symbol, ndrrOptions::list, ndrrTO::nonnegint, SeeAll::boolean, Wratio::list(nonnegint), $) :: list(list):
local i, Poss, nPoss, PPs:
description "Heuristics_EC: Selects a formulation for ECCAD using heuristics.",
			"Input: Two lists of polynomials (the first assumed equational constraints) and a variable ordering.  Optionally a heuristic choice (S, N, NS or SN) for which the default is NS.",
			"Ouput: A single ECCAD formulation or list of the ECCAD formulations suggested by the heuristic.":

Poss := Formulations_ECCAD( [ecs, necs] ):
nPoss := nops(Poss):
userinfo(2, {'ProjectionCAD'}, nPoss, "possible formulations have been identified"):
PPs := table():
for i from 1 to nPoss do 
	PPs[i] := TTI_ECPP( Poss[i], vars ):
od:
Heuristics_RunOnPPs( Poss, PPs, vars, heuristic, ndrrOptions, ndrrTO, SeeAll, Wratio, false):
end proc:

#####################################################################

Heuristics_TTICADQFF := proc( IN::list(list(polynom)), vars::list(symbol), heuristic::symbol, ndrrOptions::list, ndrrTO::nonnegint, SeeAll::boolean, splitting::boolean, Wratio::list(nonnegint), $) :: list(list):
local F, G, i, Poss, nPoss, PPs:
description "Heuristics_TTICADQFF : Selects formulation(s) for a TTICADQFF using heuristics.",
			"Input: Two lists of polynomials (the first assumed equational constraints) and a variable ordering.  Optionally a heuristic choice (S, N, NS or SN) for which the default is NS and the choice to consider splittings (for which the default is true).",
			"Ouput: A single TTICAD formulation or list of the TTICAD formulations suggested by the heuristicselected partitions.":

F := IN[1]:
G := IN[2]:
if splitting=false then
	Poss := ( map(X->[X],Formulations_ECCAD([F,G])) ):
else 
	Poss := ( Formulations_QFF( [expand(F), expand(G)] ) ):
fi:
nPoss := nops(Poss):
userinfo(2, {'ProjectionCAD'}, nPoss, "possible formulations have been identified"):
PPs := table():
for i from 1 to nPoss do 
	PPs[i] := TTI_TTIPP( Poss[i], vars ):
od:
Heuristics_RunOnPPs( Poss, PPs, vars, heuristic, ndrrOptions, ndrrTO, SeeAll, Wratio, false):
end proc:

#####################################################################

Heuristics_TTICAD := proc(LL::list(list(list(polynom))), vars::list(symbol), modular::boolean, heuristic::symbol, ndrrOptions::list, ndrrTO::nonnegint, SeeAll::boolean, splitting::boolean, Wratio::list(nonnegint), $) :: list(list):
local i, Poss, nPoss, PPs, F, G, Suggestion:
description "Heuristics_TTICAD: Selects formulation(s) for a sequence of QFFs for TTICAD using heuristics.",
			"Input: A list of lists each representing a QFF.  Each QFF has two lists of polynomials, the first assumed the equational constraints.",
					"Optionally a heuristic choice (S, N, NS or SN) for which the default is NS and the choice to consider splittings (for which the default is true).",
					"Note: it is assumed that the QFFs come from a disjunctive normal form.  I.e. trivial merging is not possible.",
			"Ouput: A single TTICAD formulation or list of the TTICAD formulations suggested by the heuristicselected partitions.":

if modular=false then
	Poss := Formulations_TTICAD(LL, splitting):
	nPoss := nops(Poss):
	userinfo(2, {'ProjectionCAD'}, nPoss, "possible formulations have been identified"):
	PPs := table():
	for i from 1 to nPoss do 
		PPs[i] := TTI_TTIPP( Poss[i], vars ):
	od:
	Heuristics_RunOnPPs( Poss, PPs, vars, heuristic, ndrrOptions, ndrrTO, SeeAll, Wratio, false):
else
	userinfo(2, {'ProjectionCAD'}, "the modular algorithm selects the formulation of each QFF in turn."):
	Suggestion:=table():
	for i from 1 to nops(LL) do
		userinfo(2, {'ProjectionCAD'}, "now considering QFF", i):
		F:=LL[i][1]:
		G:=LL[i][2]:
		Suggestion[i] := op(Heuristics_TTICADQFF( [F, G], vars, heuristic, ndrrOptions, ndrrTO, FAIL, splitting, Wratio)):
	od:
	return( [entries( Suggestion, 'nolist')] ):
fi:
end proc:

#####################################################################

Heuristics_VariableOrdering := proc( VarsList::list, IN, algorithm::symbol, greedy::boolean, heuristic::symbol, method::symbol, ndrrOptions::list, nTimeO::nonnegint, SeeAll, Wratio::list, $ ) :: list:
local PossVarOrds, i, k, t, nPoss, PPs, E, A, B, C, F, mvar, lvars, PS, Poss, selected, FirstBlock, allvars, newFirstBlock, newVarsList, BigC:
description "Heuristics_VariableOrdering: Pick a variable ordering based on a heuristic.",
			"Input: The input for a CAD algorithm, a variable list (possibly in blocks) and an algorithm selection.",
			"Ouput: Suggested variable ordering(s).":

if heuristic = 'BrownFull' then
	if VarsList::'list'(symbol) then
		return( Heuristics_Brown( VarsList, {}, IN, false, false) ):
	else
		return( Heuristics_BrownBlocks( VarsList, IN, false, false ) ):
	fi:
elif heuristic = 'BrownBasic' then
	if VarsList::'list'(symbol) then
		return( Heuristics_Brown( VarsList, {}, IN, true, false) ):
	else
		return( Heuristics_BrownBlocks( VarsList, IN, true, false ) ):
	fi:
fi:
if greedy=false then
	if VarsList::'list'(symbol) then
		PossVarOrds := Formulation_VarOrdFree(VarsList):
	elif VarsList::'list'('list'('symbol')) then
		PossVarOrds := Formulation_VarOrdRestricted(VarsList):
	fi:
	nPoss := nops(PossVarOrds):
	userinfo(2, {'ProjectionCAD'}, nPoss, "possible formulations have been identified"):
	PPs := table(): 
	if algorithm='algorithm_CADFull' then
		if method='McCallum' then 
			for i from 1 to nPoss do 
				PPs[i] := PCAD_McCProjPolys( IN, PossVarOrds[i] ):
			od:
		elif method='Collins' then 
			for i from 1 to nPoss do 
				PPs[i] := PCAD_CCADProjPolys( IN, PossVarOrds[i] ):
			od:
		else
			error("this is an unexpected error - please report to M.England@bath.ac.uk"):
		fi:
	elif algorithm='algorithm_ECCAD' then
		for i from 1 to nPoss do 
			PPs[i] := TTI_ECPP( IN, PossVarOrds[i] ):
		od:
	elif algorithm='algorithm_TTICAD' then
		for i from 1 to nPoss do 
			PPs[i] := TTI_TTIPP( IN, PossVarOrds[i] ):
		od:
	else
		error("this is an unexpected error - please report to M.England@bath.ac.uk"):
	fi:
	return( Heuristics_RunOnPPs( PossVarOrds, PPs, [], heuristic, ndrrOptions, nTimeO, SeeAll, Wratio, true) ):
else
	if algorithm='algorithm_CADFull' then
		if VarsList::'list'(symbol) then
			return( Heuristics_GVOCADFull( VarsList, {}, IN, heuristic, ndrrOptions, Wratio, false) ):
		else
			return( Heuristics_GVOCADFullBlocks( VarsList, IN, heuristic, ndrrOptions, nTimeO, Wratio, false) ):
		fi:
	elif algorithm='algorithm_ECCAD' then
		E := expand({IN[1]}):
		A := expand({IN[1], op(IN[2])}):
		if VarsList::'list'(symbol) then
			PPs := table():
			for i from 1 to nops(VarsList) do
				mvar:=VarsList[i]:
				lvars := remove(X->X=mvar, VarsList): 
				F := PCAD_PrimSet(E, mvar):
				F := PCAD_SFBasis(F, mvar):
				if nops(VarsList)=1 then
					PPs[i]:=E:
				else
					B := PCAD_PrimSet(A, mvar):
					C := PCAD_ContSet(A, mvar):
					PPs[i] := PCAD_SetFactors( TTI_ECReducedProjectionOperator( F, B, mvar, lvars) union C ):
				fi:
			od:
			Poss:=[seq( [PPs[i], VarsList[i]], i=1..nops(VarsList))]:
			selected:=Heuristics_RunOnPPs( Poss, PPs, [], heuristic, ndrrOptions, nTimeO, FAIL, Wratio, true):
			PS:=selected[1]:
			mvar := selected[2]:
			lvars := remove(X->X=mvar, VarsList): 
			lvars := Heuristics_GVOCADFull( lvars, {}, PS, heuristic, ndrrOptions, nTimeO, Wratio, false):
			return([mvar,op(lvars)]):
		else 
			FirstBlock := VarsList[1]:
			allvars:=map(X->op(X),VarsList):
			if nops(FirstBlock)=1 then
				mvar:=op(FirstBlock):
				lvars := remove(X->X=mvar, allvars): 
				F := PCAD_PrimSet(E, mvar):
				F := PCAD_SFBasis(F, mvar):
				B := PCAD_PrimSet(A, mvar):
				C := PCAD_ContSet(A, mvar):
				PS := PCAD_SetFactors( TTI_ECReducedProjectionOperator( F, B, mvar, lvars) union C ):
				lvars := Heuristics_GVOCADFullBlocks( VarsList[2..nops(VarsList)], PS, heuristic, ndrrOptions, nTimeO, Wratio, false):
				return([mvar,op(lvars)]):
			else
				PPs := table():
				for i from 1 to nops(FirstBlock) do
					mvar:=FirstBlock[i]:
					lvars := remove(X->X=mvar, allvars): 
					F := PCAD_PrimSet(E, mvar):
					F := PCAD_SFBasis(F, mvar):
					B := PCAD_PrimSet(A, mvar):
					C := PCAD_ContSet(A, mvar):
					PPs[i] := PCAD_SetFactors( TTI_ECReducedProjectionOperator( F, B, mvar, lvars) union C ):
				od:
				Poss:=[seq( [PPs[i], FirstBlock[i]], i=1..nops(FirstBlock))]:
				selected:=Heuristics_RunOnPPs( Poss, PPs, [], heuristic, ndrrOptions, nTimeO, FAIL, Wratio, true):
				PS:=selected[1]:
				mvar := selected[2]:
				newFirstBlock := remove(X->X=mvar, FirstBlock):
				newVarsList := [newFirstBlock, op(VarsList[2..nops(VarsList)])]:
				lvars := Heuristics_GVOCADFullBlocks( newVarsList, PS, heuristic, ndrrOptions, nTimeO, Wratio, false):
				return([mvar,op(lvars)]):
			fi:
		fi:
	elif algorithm='algorithm_TTICAD' then
		t := nops(IN):
		A := table():
		E := table():
		for k from 1 to t do 
			E[k] := expand( {IN[k][1]} ):
			A[k] := expand( {IN[k][1], op(IN[k][2])} ):
		od:
		if VarsList::'list'(symbol) then
			PPs := table():
			for i from 1 to nops(VarsList) do
				mvar:=VarsList[i]:
				lvars := remove(X->X=mvar, VarsList): 
				B := table():
				C := table():
				for k from 1 to t do 
					B[k] := PCAD_PrimSet(A[k], mvar):
					C[k] := PCAD_ContSet(A[k], mvar):
					F[k] := PCAD_PrimSet(E[k], mvar):
					F[k] := PCAD_SFBasis(F[k], mvar):
				od:
				BigC := `union`(seq(C[k], k=1..t)):
				PPs[i] := PCAD_SetFactors( TTI_TTIProjectionOperator(F, B, t, mvar, lvars) union BigC ):
			od:
			Poss:=[seq( [PPs[i], VarsList[i]], i=1..nops(VarsList))]:
			selected:=Heuristics_RunOnPPs( Poss, PPs, [], heuristic, ndrrOptions, nTimeO, FAIL, Wratio, true):
			PS:=selected[1]:
			mvar := selected[2]:
			lvars := remove(X->X=mvar, VarsList): 
			lvars := Heuristics_GVOCADFull( lvars, {}, PS, heuristic, ndrrOptions, nTimeO, Wratio, false):
			return([mvar,op(lvars)]):
		else 
			FirstBlock := VarsList[1]:
			allvars:=map(X->op(X),VarsList):
			if nops(FirstBlock)=1 then
				mvar:=op(FirstBlock):
				lvars := remove(X->X=mvar, allvars):
				B := table():
				C := table():
				for k from 1 to t do 
					B[k] := PCAD_PrimSet(A[k], mvar):
					C[k] := PCAD_ContSet(A[k], mvar):
					F[k] := PCAD_PrimSet(E[k], mvar):
					F[k] := PCAD_SFBasis(F[k], mvar):
				od:
				BigC := `union`(seq(C[k], k=1..t)):
				PS := PCAD_SetFactors( TTI_TTIProjectionOperator(F, B, t, mvar, lvars) union BigC ):
				lvars := Heuristics_GVOCADFullBlocks( VarsList[2..nops(VarsList)], PS, heuristic, ndrrOptions, nTimeO, Wratio, false):
				return([mvar,op(lvars)]):
			else
				PPs := table():
				for i from 1 to nops(FirstBlock) do
					mvar:=FirstBlock[i]:
					lvars := remove(X->X=mvar, allvars): 
				B := table():
				C := table():
				for k from 1 to t do 
					B[k] := PCAD_PrimSet(A[k], mvar):
					C[k] := PCAD_ContSet(A[k], mvar):
					F[k] := PCAD_PrimSet(E[k], mvar):
					F[k] := PCAD_SFBasis(F[k], mvar):
				od:
				BigC := `union`(seq(C[k], k=1..t)):
				PPs[i] := PCAD_SetFactors( TTI_TTIProjectionOperator(F, B, t, mvar, lvars) union BigC ):
				od:
				Poss:=[seq( [PPs[i], FirstBlock[i]], i=1..nops(FirstBlock))]:
				selected:=Heuristics_RunOnPPs( Poss, PPs, [], heuristic, ndrrOptions, nTimeO, FAIL, Wratio, true):
				PS:=selected[1]:
				mvar := selected[2]:
				newFirstBlock := remove(X->X=mvar, FirstBlock):
				newVarsList := [newFirstBlock, op(VarsList[2..nops(VarsList)])]:
				lvars := Heuristics_GVOCADFullBlocks( newVarsList, PS, heuristic, ndrrOptions, nTimeO, Wratio, false):
				return([mvar,op(lvars)]):
			fi:
		fi:
	else
		error("this is an unexpected error - please report to M.England@bath.ac.uk"):
	fi:
fi:
end proc: #Heuristics_VariableOrdering

#####################################################################

Heuristics_GVOCADFull := proc( vars::list(symbol), lvars::set(symbol), P::set(polynom), heuristic::symbol, ndrrOptions::list, nTimeO::nonnegint, Wratio::list(nonnegint), retpols::boolean, $ ) :: list:
local i, cont, othervars, Avars, UAvars, PPs, PS, var, pset, Poss, selected:
description "Heuristics_GVOBlock: Use the greedy method to choose a variable ordering (based on CADFull-McCallum).",
			"Input: Variables, polynomials, and the usual optional heuristic arguments.",
			"Ouput: Suggested variable ordering.":

PS:=P:
Avars:=[]:
UAvars:=vars:
while nops(UAvars)>1 do
	PPs := table():
	for i from 1 to nops(UAvars) do
		var:=UAvars[i]:
		othervars := convert(remove(X->X=var,UAvars), 'set') union lvars:
		cont:=PCAD_ContSet(PS,var):  
		pset:=PCAD_McCallumProj( PS, var, convert(othervars, 'list') ):
		PPs[i]:=PCAD_SetFactors( pset union cont ):
	od:
	Poss:=[seq( [PPs[i], UAvars[i]], i=1..nops(UAvars))]:
	selected:=Heuristics_RunOnPPs( Poss, PPs, [], heuristic, ndrrOptions, nTimeO, FAIL, Wratio, true):
	Avars:=[op(Avars), selected[2]]:
	UAvars:=remove(X->X=selected[2], UAvars):
	PS:=selected[1]:
od:
if nops(UAvars)<>1 then 
	error("This is an unexpected error.  Please report to M.England@bath.ac.uk"):
else
	var := op(UAvars):
	cont := PCAD_ContSet(PS,var):  
	pset := PCAD_McCallumProj( PS, var, convert(lvars, 'list') ):
	PS := PCAD_SetFactors( pset union cont ):
	Avars:=[op(Avars), op(UAvars)]:
fi:
if retpols=false then
	return(Avars):
else
	return(Avars,PS):
fi:
end proc:

#####################################################################

Heuristics_GVOCADFullBlocks := proc( VarsBlocks::list(list(symbol)), P::set(polynom), heuristic::symbol, ndrrOptions::list, nTimeO::nonnegint, Wratio::list(nonnegint), retpols::boolean, $ ) :: list:
local i, k, numBlocks, PS, Avars, lvars, vars:
description "Heuristics_GVOBlock: Use the greedy method to choose a variable ordering (based on CADFull-McCallum) when variables in blocks.",
			"Input: Variables, polynomials, and the usual optional heuristic arguments.",
			"Ouput: Suggested variable ordering.":

numBlocks:=nops(VarsBlocks):
PS := P:
Avars := []:
for k from 1 to numBlocks do
	lvars:={seq(op(VarsBlocks[i]), i=(k+1)..numBlocks)}:
	vars,PS := Heuristics_GVOCADFull( VarsBlocks[k], lvars, PS, heuristic, ndrrOptions, nTimeO, Wratio, true):
	Avars:=[op(Avars), op(vars)]:
od:
if retpols=false then
	return(Avars):
else
	return(Avars,PS):
fi:
end proc:

#####################################################################

Heuristics_Brown := proc( VarsList::list(symbol), lvars::set(symbol), P::set(polynom), basic::boolean, retpols::boolean, $ ) :: list:
local i, nPoss, var, UAvars, Avars, BH1s, BH2s, BH3s, PS, BHF, Poss, Suggestions, SuggInd:
description "Heuristics_Brown: Pick a variable ordering for CADFull based on Brown's heuristic.  This is coded separetly as unlike other heuristics it does not calculate projection polynomials before making the variable choice.",
			"Input: The list of variables, the list of input polynomials to CADFull.  If basic then only the input is studied while if full the projection layers are studied in turn.",
			"Ouput: Suggested variable ordering(s).":

PS:=P:
Avars:=[]:
UAvars:=VarsList:
while nops(UAvars)>1 do
	Poss := UAvars:
	nPoss := nops(Poss):
	BHF := table():
	for i from 1 to nPoss do
		var:=UAvars[i]:
		BHF[i] := Heuristics_BHFeatures(PS, var);
	od:
	BH1s := [seq( BHF[i][1], i = 1..nPoss)]: 
	BH2s := [seq( BHF[i][2], i = 1..nPoss)]: 
	BH3s := [seq( BHF[i][3], i = 1..nPoss)]: 
	Suggestions := Heuristics_PickSuggestions( Poss, BH1s, 'pick'='maximise'):
	if nops(Suggestions) = 1 then
		var := op(Suggestions):
	else
		SuggInd := select( X -> has( map(Y->is(Y=Poss[X]),Suggestions) ,true)=true, [seq(1..nops(Poss))] ):
		Poss := Poss[SuggInd]:
		BH2s := BH2s[SuggInd]:
		BH3s := BH3s[SuggInd]:
		Suggestions := Heuristics_PickSuggestions( Poss, BH2s, 'pick'='maximise'):
		if nops(Suggestions) = 1 then
			var := op(Suggestions):
		else
			SuggInd := select( X -> has( map(Y->is(Y=Poss[X]),Suggestions) ,true)=true, [seq(1..nops(Poss))] ):
			Poss := Poss[SuggInd]:
			BH3s := BH3s[SuggInd]:
			Suggestions := Heuristics_PickSuggestions( Poss, BH3s, 'pick'='maximise'):
			if nops(Suggestions) = 1 then
				var := op(Suggestions):
			else
				var := Suggestions[1]:
			fi:
		fi:
	fi:
	Avars := [op(Avars), var ]:
	UAvars := remove(X->X=var, UAvars):
	if basic=false then
		PS := PCAD_McCallumProj( PS, var, [op(lvars), op(UAvars)] ):
	else
		PS := PS:
	fi:
od:
if nops(UAvars)<>1 then 
	error("This is an unexpected error.  Please report to M.England@bath.ac.uk"):
else
	var := op(UAvars):
	Avars:=[op(Avars), op(UAvars)]:
	if basic=false then
		PS := PCAD_McCallumProj( PS, var, convert(UAvars, 'list') ):
	else
		PS := PS:
	fi:
fi:
if retpols=false then
	return(Avars):
else
	return(Avars,PS):
fi:
end proc:

#####################################################################

Heuristics_BrownBlocks := proc( VarsBlocks::list(list(symbol)), P::set(polynom), basic::boolean, retpols::boolean, $ ) :: list:
local i, k, numBlocks, PS, Avars, lvars, vars:
description "Heuristics_BrownBlock: Pick a variable ordering for CADFull based on Brown's heuristic.  The case when the variables are in blocks.",  
			"Input: The list of variable blocks, the list of input polynomials to CADFull.",
			"Ouput: Suggested variable ordering(s).":

numBlocks := nops(VarsBlocks):
PS := P:
Avars := []:
for k from 1 to numBlocks do
	lvars:={seq(op(VarsBlocks[i]), i=(k+1)..numBlocks)}:
	vars,PS := Heuristics_Brown( VarsBlocks[k], lvars, PS, basic, true):
	Avars:=[op(Avars), op(vars)]:
od:
if retpols=false then
	return(Avars):
else
	return(Avars,PS):
fi:
end proc:

#####################################################################
###  Section 9: SubCAD (David Wilson - D.J.Wilson@bath.ac.uk
#####################################################################

LCAD_ProjCADLiftOneLayered:=proc( ProjPolys::set, vars::list, meth, finalCAD, retcad, out, failure, $ )
    local i, j, k, switch, dim, R, contset, pset, cad, tmpcad, tmppset, cell, SPrc, SPcube, cellIndex, alpha, falpha, DP, ICZD:
    description "LCAD_ProjCADLiftOneLayered: Constructs a CAD lifting only full dimensional cells to produce a 1-layered CAD.":

    if kernelopts('opaquemodules')=true then
        switch:=0:
        kernelopts('opaquemodules'=false):
    fi:
    dim:=nops(vars):
    contset:={}:
    for i from 1 to dim do
        R[i]:=RegularChains:-PolynomialRing( vars[-i..-1] ):
        pset[dim+1-i]:=select(has, remove(has,ProjPolys,vars[1..i-1]), vars[i]):
        pset[dim+1-i]:=pset[dim+1-i] union contset:
        contset:=PCAD_ContSet(pset[dim+1-i], vars[i]):
        pset[dim+1-i]:=PCAD_PrimSet(pset[dim+1-i], vars[i]):
        pset[dim+1-i]:=convert(pset[dim+1-i],'list'):
        pset[dim+1-i]:=convert(PCAD_SetFactors(pset[dim+1-i]),'list'):
        if pset[i]=[] then WARNING("there are no projection polynomials in %1", R[i]['variables']): fi:
    od:
    pset[1]:=remove(X -> sturm(X, vars[-1], -infinity, infinity) = 0, pset[1]):
    cad[1]:=CADFull(pset[1], vars[-1..-1],'method'=meth,'output'=out):
#    cad[1]:=map(X->X[1], cad[1]):
    userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"produced CAD of", R[1]['variables'], "-space with", nops(cad[1]), "cells"):
    if retcad=1 then
        return(cad[1]):
    fi:
    for i from 2 to dim do
        tmpcad:=table():
        for j from 1 to nops(cad[i-1]) do
            cell:=cad[i-1][j]:
            cellIndex:=cell[1]:
            #Check the cells for full dimensionality
            if PCAD_IsCellFullDim(cellIndex) then # new
                SPrc:=op(1,op(-1,cell)):
                SPcube:=op(2,op(-1,cell)):
                alpha:=PCAD_SPtoRootOf( SPrc,SPcube,R[i-1] ):
                tmppset:=[]:
                for k from 1 to nops(pset[i]) do
                    falpha:=subs(alpha,op(k,pset[i])):
                    falpha:=is(falpha=0):
                    if falpha=true then
                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"a projection polynomial was nullified on cell", cellIndex):
                        ICZD:=PCAD_IsCellZeroDim( cellIndex ):
                        if i=dim then
                            if finalCAD='SI' then
                                tmppset:=[op(tmppset)]:
                                userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"It was discarded since this is the final lift and only sign invariance is required."):
                            else
                                if ICZD=true then
                                    DP:=PCAD_MinimalDelineatingPolynomial(op(k,pset[i]),alpha,R[i]):
                                    if DP::constant then
                                        tmppset:=[op(tmppset)]:
                                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"It was discarded since cell", cellIndex, "is zero-dim and the minimal delineating polynomial is a constant."):
                                    else
                                        tmppset:=[op(tmppset), DP]:
                                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "Since cell", cellIndex, "is zero-dim a delineating polynomial was used."):
                                        userinfo(4, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "The delineating polynomial used is", DP):
                                    fi:
                                else
                                    userinfo(4, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "There is a problem.  The polynomial", op(k,pset[i]), "is nullified by", alpha, "which is a cell of dim>0."):
                                    if failure='err' then
                                        error("there is nullification on cell %1 which has dim>0. The outputted CAD cannot be guaranteed order-invariant.", cellIndex):
                                    elif failure='giveFAIL' then
                                        return(FAIL):
                                    else
                                        WARNING("there is nullification on cell %1 which has dim>0. The nullified polynomial was discarded. The outputted CAD may not be order-invariant.", cellIndex):
                                        tmppset:=[op(tmppset)]:
                                    fi:
                                fi:
                            fi:
                        else
                            if meth='McCallum' then
                                if ICZD=true then
                                    DP:=PCAD_MinimalDelineatingPolynomial(op(k,pset[i]),alpha,R[i]):
                                    if DP::constant then
                                        tmppset:=[op(tmppset)]:
                                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"It was discarded since cell", cellIndex, "is zero-dim and the minimal delineating polynomial is a constant."):
                                    else
                                        tmppset:=[op(tmppset), DP]:
                                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"Since cell", cellIndex, "is zero-dim a delineating polynomial was used."):
                                        userinfo(4,'PCAD_ProjCADLift',"The delineating polynomial used is", DP):
                                    fi:
                                else
                                    userinfo(4, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "There is a problem.  The polynomial", op(k,pset[i]), "is nullified by", alpha, "which is a cell of dim>0."):
                                    if failure='err' then
                                        error("there is nullification on cell %1 which has dim>0. The outputted CAD cannot be guaranteed sign-invariant.", cellIndex):
                                    elif failure='giveFAIL' then
                                        return(FAIL):
                                    else
                                        WARNING("there is nullification on cell %1 which has dim>0. The nullified polynomial was discarded. The outputted CAD may not be sign-invariant.", cellIndex):
                                        tmppset:=[op(tmppset)]:
                                    fi:
                                fi:
                            else
                                tmppset:=[op(tmppset)]:
                                userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"The nullified polynomial on cell", cellIndex, "was discarded since we are using Collins lifting."):
                            fi:
                        fi:
                    else tmppset:=[op(tmppset),op(k,pset[i])]:
                    fi:
                od:
                #tmpcad[j]:=op(PCAD_GenerateStack(
                # cell, tmppset, R[i], out )): # Incorrect method of lifting
                tmpcad[j]:=op(CADGenerateStack( cell, tmppset, vars[-i..-1],'output'= out )):

                else
                    #If not satisfying the dimension requirements, replace
                    # by something trivial that fails the test
                    tmpcad[j]:=[[2],[]]:
                :fi #new
           od:
        cad[i]:=[seq( tmpcad[k], k=1..nops(cad[i-1]) )]:
        #cad[i]:=map(X->op(1,X),cad[i]): #Not needed with correct liting
        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"produced CAD of", R[i]['variables'], "-space with", nops(cad[i]), "cells"):
        if retcad=i then
            return(cad[i]):
        fi:
    od:
    if switch=0 then
        kernelopts('opaquemodules'=true):
    fi:
    # We still need to select only the full dimensional cells from the
    # top layer
    cad[dim]:=select(LCAD_IsCellRepOneLayeredDim,cad[dim]):
    if out='list' or out='listwithrep' then
        return( cad[dim] ):
    elif out='rootof' then
        tmpcad:=map(X -> X[2],cad[dim]):
        tmpcad:=map(Y-> map(X -> `if`(op(0,X)=And, op(X), X), Y), tmpcad);
        return(tmpcad):
    elif out='piecewise' then
        return( PCAD_LWRCADtoPWCAD( cad[dim] ) ):
    fi:
end proc:

#####################################################################

LCAD_ProjCADLiftNLayered:=proc( ProjPolys::set,n, vars::list, meth, finalCAD, retcad, out, failure, $ )
    local i, j, k, switch, dim, R, contset, pset, cad, tmpcad, tmppset, cell, SPrc, SPcube, cellIndex, alpha, falpha, DP, ICZD:
    description "PCAD_ProjCADLiftNLayered: Given a list of projection polynomials and a list of variables, compute an n-layered subCAD.":

    if kernelopts('opaquemodules')=true then
        switch:=0:
        kernelopts('opaquemodules'=false):
    fi:
    dim:=nops(vars):
    contset:={}:
    for i from 1 to dim do
        R[i]:=RegularChains:-PolynomialRing( vars[-i..-1] ):
        pset[dim+1-i]:=select(has, remove(has,ProjPolys,vars[1..i-1]), vars[i]):
        pset[dim+1-i]:=pset[dim+1-i] union contset:
        contset:=PCAD_ContSet(pset[dim+1-i], vars[i]):
        pset[dim+1-i]:=PCAD_PrimSet(pset[dim+1-i], vars[i]):
        pset[dim+1-i]:=convert(pset[dim+1-i],'list'):
        pset[dim+1-i]:=convert(PCAD_SetFactors(pset[dim+1-i]),'list'):
        if pset[i]=[] then WARNING("there are no projection polynomials in %1", R[i]['variables']): fi:
    od:
    pset[1]:=remove(X -> sturm(X, vars[-1], -infinity, infinity) = 0, pset[1]):
    cad[1]:=CADFull( pset[1], vars[-1..-1], 'method'=meth, 'output'=out ):
    userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"produced CAD of", R[1]['variables'], "-space with", nops(cad[1]), "cells"):
    if retcad=1 then
        return(cad[1]):
    fi:
    for i from 2 to dim do
        tmpcad:=table():
        for j from 1 to nops(cad[i-1]) do
            cell:=cad[i-1][j]:
            cellIndex:=cell[1]:
            #Check the dimension requirements
            if LCAD_IsCellNLayeredDim(cellIndex,n) then # new
                SPrc:=op(1,op(-1,cell)):
                SPcube:=op(2,op(-1,cell)):
                alpha:=PCAD_SPtoRootOf( SPrc,SPcube,R[i-1] ):
                tmppset:=[]:
                for k from 1 to nops(pset[i]) do
                    falpha:=subs(alpha,op(k,pset[i])):
                    falpha:=is(falpha=0):
                    if falpha=true then
                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"a projection polynomial was nullified on cell", cellIndex):
                        ICZD:=PCAD_IsCellZeroDim( cellIndex ):
                        if i=dim then
                            if finalCAD='SI' then
                                tmppset:=[op(tmppset)]:
                                userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"It was discarded since this is the final lift and only sign invariance is required."):
                            else
                                if ICZD=true then
                                    DP:=PCAD_MinimalDelineatingPolynomial(op(k,pset[i]),alpha,R[i]):
                                    if DP::constant then
                                        tmppset:=[op(tmppset)]:
                                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"It was discarded since cell", cellIndex, "is zero-dim and the minimal delineating polynomial is a constant."):
                                    else
                                        tmppset:=[op(tmppset), DP]:
                                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "Since cell", cellIndex, "is zero-dim a delineating polynomial was used."):
                                        userinfo(4, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "The delineating polynomial used is", DP):
                                    fi:
                                else
                                    userinfo(4, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "There is a problem.  The polynomial", op(k,pset[i]), "is nullified by", alpha, "which is a cell of dim>0."):
                                    if failure='err' then
                                        error("there is nullification on cell %1 which has dim>0. The outputted CAD cannot be guaranteed order-invariant.", cellIndex):
                                    elif failure='giveFAIL' then
                                        return(FAIL):
                                    else
                                        WARNING("there is nullification on cell %1 which has dim>0. The nullified polynomial was discarded. The outputted CAD may not be order-invariant.", cellIndex):
                                        tmppset:=[op(tmppset)]:
                                    fi:
                                fi:
                            fi:
                        else
                            if meth='McCallum' then
                                if ICZD=true then
                                    DP:=PCAD_MinimalDelineatingPolynomial(op(k,pset[i]),alpha,R[i]):
                                    if DP::constant then
                                        tmppset:=[op(tmppset)]:
                                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"It was discarded since cell", cellIndex, "is zero-dim and the minimal delineating polynomial is a constant."):
                                    else
                                        tmppset:=[op(tmppset), DP]:
                                        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"Since cell", cellIndex, "is zero-dim a delineating polynomial was used."):
                                        userinfo(4,'PCAD_ProjCADLift',"The delineating polynomial used is", DP):
                                    fi:
                                else
                                    userinfo(4, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'}, "There is a problem.  The polynomial", op(k,pset[i]), "is nullified by", alpha, "which is a cell of dim>0."):
                                    if failure='err' then
                                        error("there is nullification on cell %1 which has dim>0. The outputted CAD cannot be guaranteed sign-invariant.", cellIndex):
                                    elif failure='giveFAIL' then
                                        return(FAIL):
                                    else
                                        WARNING("there is nullification on cell %1 which has dim>0. The nullified polynomial was discarded. The outputted CAD may not be sign-invariant.", cellIndex):
                                        tmppset:=[op(tmppset)]:
                                    fi:
                                fi:
                            else
                                tmppset:=[op(tmppset)]:
                                userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"The nullified polynomial on cell", cellIndex, "was discarded since we are using Collins lifting."):
                            fi:
                        fi:
                    else tmppset:=[op(tmppset),op(k,pset[i])]:
                    fi:
                od:
                #tmpcad[j]:=op(PCAD_GenerateStack(
                # cell, tmppset, R[i], out )): # incorrect method of
                # lifting
                tmpcad[j]:=op(CADGenerateStack( cell, tmppset, vars[-i..-1], 'output'=out )):
                else
                    #If sub-dimensional, add a cell that will never
                    # satisfy the lifting property
                    tmpcad[j]:=[[2$n],[]]:
                :fi #new
           od:
        cad[i]:=[seq( tmpcad[k], k=1..nops(cad[i-1]) )]:
        #cad[i]:=map(X->op(1,X),cad[i]): #Not needed now using correct lfiting
        userinfo(3, {'PCAD_ProjCADLift', 'CADLifting', 'CADFull'},"produced CAD of", R[i]['variables'], "-space with", nops(cad[i]), "cells"):
        if retcad=i then
            return(cad[i]):
        fi:
    od:
    if switch=0 then
        kernelopts('opaquemodules'=true):
    fi:
    #New code to remove top level cells not of the right dimension
    tmpcad:=cad[dim]:
    cad[dim]:=[]:
    #Cycle through and remove as needed
    for i from 1 to nops(tmpcad) do
        if tmpcad[i]=[[2$n],[]] then
        else
            if LCAD_IsCellNLayeredDim(tmpcad[i][1],n) then
                cad[dim]:=[op(cad[dim]),tmpcad[i]]:
            fi:
        fi:
    od:
    if out='list' or out='listwithrep' then
        return( cad[dim] ):
    elif out='rootof' then
        tmpcad:=map(X -> X[2],cad[dim]):
        tmpcad:=map(Y-> map(X -> `if`(op(0,X)=And, op(X), X), Y), tmpcad);
        return(tmpcad):
    elif out='piecewise' then
        return( PCAD_LWRCADtoPWCAD( cad[dim] ) ):
    fi:
end proc:

#####################################################################

LCAD_IsCellNLayeredDim:=proc( cellindex::list, n, $ ) :: boolean:
    local modtot,N,i:
    description "PCAD_IsCellFullDimPlusOne: Check to see if a cell index indicates it is of sufficient dimension to contribute to the top n layers of a CAD.":
    N:=nops(cellindex):
    #Calculate dimension by summing parities of indices
    modtot:=add(i mod 2, i in cellindex):
    #Top n layers have dimension >= N+1-n
    if modtot >= (N+1-n) then
        return(true):
    elif modtot<(N+1-n) then
        return(false):
    else
        error("could not determine parity. This is an unexpected error message - please report to D.J.Wilson@bath.ac.uk"):
    fi:
end proc:

#####################################################################

LCAD_IsCellRepOneLayeredDim:=proc(cell) :: boolean:
    description "LCAD_IsCellRepOneLayeredDim: returns true if cell is full dimensional and false otherwise":
      if LCAD_IsCellNLayeredDim(cell[1],1) then
        return true:
      else
        return false:
      fi:
end proc:

#####################################################################

LCAD_IsCellRepNLayeredDim:=proc(cell,n) :: boolean:
    description "LCAD_IsCellRepNLayeredDim: Check to see if a cell is in the top n layers of its CAD.":
    #Calculate by a call to LCAD_IsCellNLayeredDim
    if LCAD_IsCellNLayeredDim(cell[1],n) then
        return true:
    else
        return false:
    fi:
end proc:

#####################################################################

LCAD_GenerateSeparateProj:=proc(F,vars, {method::symbol:='McCallum'})
    local n,Proj,ProjAll,i:
    description "LCAD_GenerateSeparateProj: Create a list of projection polynomials sorted according to their main variable.":

    n:=nops(vars):
    Proj:=[{}$n]:
    if method='McCallum' then
        ProjAll:=(PCAD_McCProjPolys(F,vars)):
    elif method='Collins' then
        ProjAll:=(PCAD_CCADProjPolys(F,vars)):
    else
        print("Method not specified or unrecognised. McCallum's projection operator will be used"):
        ProjAll:=(PCAD_McCProjPolys(F,vars)):
    fi:
    for i from 1 to n-1 do
      Proj[i]:=convert(select(a->has(a,vars[n-i+1]),remove(a->has(a,vars[n-i]),ProjAll)),list):
    end do:
    Proj[n]:=convert(select(a->has(a,vars[1]),ProjAll),list):
    return(Proj):
end proc:

#####################################################################

LCAD_LWRCADtoPWCAD:=proc( LCAD::list ) :: 'piecewise':
    local dim, Tree, Leaf, i, NumRoots, SubCAD, Rep, j:
    description "LCAD_LWRCADtoPWCAD: Convert a layered CAD in listwithrep output format to a piecewise layeredCAD.":

    dim:=convert(map(X->nops(X[1]),LCAD),set):
    if nops(dim)<>1 then
        error("cells have indices of different length. This is an unexpected error message - please report to D.J.Wilson@bath.ac.uk"):
    else
        dim:=op(dim):
    fi:
    if dim=0 then
        if nops(LCAD)<>1 or LCAD[1][1]<>[] or LCAD[1][2]<>[] then
            error("converting to piecewise did not go as expected. This is an unexpected error message - please report to D.J.Wilson@bath.ac.uk"):
        fi:
        Leaf:=LCAD[1][3]:
        return(Leaf):
    fi:
    NumRoots:=max(map(X->X[1,1],LCAD)):
    Tree:=[]:
    for i from 1 to NumRoots do
        SubCAD:=select(X->X[1][1]=i, LCAD):
        Rep:=convert(map(X->X[2,1],SubCAD),set):
        if nops(Rep)>1 then
            error("cells on same branch have different representations. This is an unexpected error message - please report to D.J.Wilson@bath.ac.uk"):
        elif nops(Rep)=0 then
            Rep:='branch'='truncated':
            SubCAD:=[op(SubCAD), [[i,seq(1+j-j,j=1..dim-1)], [seq(Rep+(j-j),j=1..dim)], ["**************"]] ]:
        else
            Rep:=op(Rep):

        fi:
            SubCAD:=map(X->[X[1][2..dim], X[2][2..dim], X[3]], SubCAD):
            Tree:=[op(Tree), Rep, LCAD_LWRCADtoPWCAD( SubCAD )]:
    od:
    if nops(Tree)=2 then
        return(Tree[2])
    else
        piecewise(op(Tree)):
    fi:
end proc:

#####################################################################

LCAD_Distribution:=proc( F::{list,set}, vars::list, method::symbol ) :: 'list':
    local pset,FC,i,DL,cind,cdim,n:
    description "LCAD_Distribution: Generates a list with the number of cells of each dimension in a full CAD using McCallum projection.",
                "Input: A list or set of polynomials and a list of variables.",
                "Output: A list of the number of cells of each dimension.":
    #Initiate a list to store the number of cells of each dimension
    DL:=[0$(nops(vars)+1)]:
    if method='McCallum' then 
        pset := PCAD_McCProjPolys(convert(F,set),vars):
        userinfo(2, {'ProjectionCAD', 'CADFull'}, "produced set of", nops(pset), "projection factors using the", method, "algorithm."):
        userinfo(4, {'ProjectionCAD', 'CADFull'}, "the projection factors are:", pset):
    elif method='Collins' then 
        pset := PCAD_CCADProjPolys(convert(F,set),vars):
        userinfo(2, {'ProjectionCAD', 'CADFull'}, "produced set of", nops(pset), "projection factors using the", method, "algorithm."):
        userinfo(4, {'ProjectionCAD', 'CADFull'}, "the projection factors are:", pset):
    fi: 
    FC:=PCAD_ProjCADLift( pset, vars, method, 'SI', 0, list, 'warn' ):
    for i from 1 to nops(FC) do
      #Extract the cell index
      cind:=FC[i][1]:
      #Calculate the dimension by summing the parities of the indices
      cdim:=add(n mod 2, n in cind);
      #Add 1 to the appropriate entry (cdim+1) in the list
      DL[cdim+1]:=DL[cdim+1]+1:
    end do:
    return(DL):
end proc:

#####################################################################

LCAD_NormDistribution:=proc( F::{list,set}, vars::list, method::symbol ) :: 'list':
    local tmpres, tot, n:
    description "LCAD_NormDistribution: Generates a list with the normed number of cells of each dimension in a full CAD using McCallum/Collins projection.",
                "Input: A list or set of polynomials and a list of variables.",
                "Output: A list of the normed number of cells of each dimension.":

    #Compute the number of cells using CADDistribution
    tmpres:=LCAD_Distribution(F,vars, method):
    #Calculate the total cells
    tot:=add(n, n in tmpres):
    #Divide each entry by the total
    return( tmpres/tot );
end proc:

#####################################################################

LCAD_Recursive:=proc(F::{list,set},vars::list,C::list,LD::list,meth::symbol,outp::symbol,resetprojection::boolean, $)
    global PolySets:
    local Fset,n,i, numlays,R,Rec,CAD1,CAD1Full,CAD1Rec,RecCAD,CurrCADFull,CurrCADRec,cell,cellstack,CADLFull,CADLRec,CADlayer,switch:
    description "LCADRecursive: Computes a recursive layered CAD in list format. ",
                "Input: A set of polynomials, list of ordered variables, list of cells a previous recursive call terminated at (possibly empty), a layered CAD in list format from a previous recursive call (possibly empty). Also, an optional parameter stating whether to reset the projection polynomials.",
                "Output: A layered CAD in list format of one more layer than inputted (or 1 layer if an empty CAD is given), an unevaluated recursive call to produce the following layer (using the 'value' call).":

    if kernelopts('opaquemodules')=true then 
        switch:=0: 
        kernelopts('opaquemodules'=false): 
    fi:
    Fset:=convert(F,set):
    #Check if projection polynomials have already been defined
    if type(PolySets,list) then
      userinfo(3, 'LCADRecursive', "Projection polynomials already defined.");
      #If the projection polys have been defined then we need to check the resetproj boolean to see if
      # they need recalculated
      #  if convert(F,set) subset `union`(seq(convert(PolySets[i],set),i=1..nops(PolySets))) then
      if not resetprojection then
        userinfo(3, 'LCADRecursive', "Projection polynomials left as defined.");
      else
        userinfo(3, 'LCADRecursive', "Projection polynomials to be reset as resetproj defined to be true.");
        PolySets:=LCAD_GenerateSeparateProj(Fset,vars,meth):
      end if:
    else
    #If none defined, recompute using LCAD_GenerateSeparateProj
      userinfo(3, 'LCADRecursive', "No projection polynomials previously defined. Recalculating.");
      PolySets:=LCAD_GenerateSeparateProj(Fset,vars,meth):
    end if:
    n:=nops(vars):
    #Calculate the number of layers required
    if nops(C)=0 then
      numlays:=1:
    else
    #Compute numlays by looking at the first cell index and counting the
    # number of even indices, and adding 1.
      numlays:=1+add((i+1) mod 2,i in C[1][1]):   
    fi:
    userinfo(3, 'LCADRecursive', "calculated that a", numlays, "layered CAD should be produced.");
    #Create polynomial rings for each layer and split the truncated cell
    # list according to the dimension of each cell
    for i from 1 to n do
      R[i]:=RegularChains:-PolynomialRing(vars[-i..-1]):
        Rec[i]:=select(a->nops(a[1])=i,C):   
    end do:
    #Treat the one layered case separately as we have to do the initial
    # splitting of R^1 and dimension checking is slightly simpler.
    if numlays=1 then
      #Generate a stack over the empty cell (i.e. split R^1)
      CAD1:=CADFull(PolySets[1],vars[-1..-1], 'method'=meth,'output'=list):
      #Separate into full dimensional and deficient cells - store in a
      # temporary variable
      CAD1Full:=select(a->LCAD_IsCellRepOneLayeredDim(a),CAD1);
      CAD1Rec:=remove(a->LCAD_IsCellRepOneLayeredDim(a),CAD1);
      RecCAD:=CAD1Rec:
      CurrCADFull:=CAD1Full:
      CurrCADRec:=CAD1Rec:
      userinfo(3, 'LCADRecursive',"produced", nops(CurrCADFull)+nops(CurrCADRec), "cells of", R[1]['variables'], "-space."):
      userinfo(3, 'LCADRecursive',"discarded", nops(CurrCADRec), "deficient cells of", R[1]['variables'], "-space."):
      for i from 2 to max(2,n) do
        CADLFull:=[]:
        CADLRec:=[]:
        for cell in CurrCADFull do
          #Take each cell in the current list of full cells. Build a stack
          # and separate into full and deficient cells.
          cellstack:=CADGenerateStack(cell,PolySets[i],vars[-i..-1],'output'=list);
          CADLFull:=[op(CADLFull),op(select(a->LCAD_IsCellRepOneLayeredDim(a),cellstack))]:
          CADLRec:=[op(CADLRec),op(remove(a->LCAD_IsCellRepOneLayeredDim(a),cellstack))]:
        end do:
        #Update temporary variables
        CurrCADFull:=CADLFull:
        CurrCADRec:=CADLRec:
        userinfo(3, 'LCADRecursive',"produced", nops(CurrCADFull)+nops(CurrCADRec), "cells of", R[i]['variables'], "-space."):
        userinfo(3, 'LCADRecursive',"discarded", nops(CurrCADRec), "deficient cells of", R[i]['variables'], "-space."):
        #Add all deficient cells to the list of truncated cells to be
        # returned in the recursive call.
        RecCAD:=[op(RecCAD),op(CurrCADRec)]:
      end do:
      #Return the layered CAD along with a recursive call (using the value
      # symbol '%') for the next layer.
      return CurrCADFull,%LCADRecursive(F,vars,RecCAD,CurrCADFull,'method'=meth,'output'=outp,'resetproj'=false):
    fi:
    #Now deal with the case where we already have a layered CAD and are
    # building over the top of it.
    #Store the previous layered CAD in the final output CAD variable (to
    # be updated later)
    CADlayer:=LD;
    #We initialise the variables to be empty as they will be populated by
    # the bottom layer truncated cell after the first loop of the
    # following for loop.
    CurrCADFull:=[]:
    RecCAD:=[]:
    for i from 1 to n do
      CADLFull:=[]:
      CADLRec:=[]:
      for cell in CurrCADFull do
        #For each cell in the current layer generate a stack using the
        # appropriate projection polynomials and sort into cells with high
        # enough dimension and those deficient.
          cellstack:=CADGenerateStack(cell,PolySets[i],vars[-i..-1],'output'=list);
        CADLFull:=[op(CADLFull),op(select(a->LCAD_IsCellRepNLayeredDim(a,numlays),cellstack))]:
        CADLRec:=[op(CADLRec),op(remove(a->LCAD_IsCellRepNLayeredDim(a,numlays),cellstack))]:
      end do:
      #Update the current CAD level with those cells generated and those
      # truncated cells inputted into the procedure of the correct dimension.
      CurrCADFull:=CADLFull:
      userinfo(3, 'LCADRecursive',"produced", nops(CurrCADFull)+nops(CurrCADRec), "cells of", R[i]['variables'], "-space."):
      CurrCADFull:=[op(CurrCADFull),op(Rec[i])]:
      CurrCADRec:=CADLRec:
      userinfo(3, 'LCADRecursive',"discarded", nops(CurrCADRec), "deficient cells of", R[i]['variables'], "-space."):
      #Add all deficient cells to the list of truncated cells to be
      # returned in the recursive call.
      RecCAD:=[op(RecCAD),op(CurrCADRec)]:
    end do:
    #Combine the layered CAD inputted into the procedure with newly
    # computed cells.
    CADlayer:=[op(CADlayer),op(CurrCADFull)]:

    if switch=0 then 
        kernelopts('opaquemodules'=true): 
    fi:
    #Return the new layered CAD along with a recursive call (using the
    # value symbol '%') for the next layer.
    return CADlayer,%LCADRecursive(F,vars,RecCAD,CADlayer,'method'=meth,'output'=outp,'resetproj'=false):
end proc:

#####################################################################

LCAD_RecursiveLR:=proc(F::{list,set},vars::list,C::list,LD::list,meth::symbol,outp::symbol,resetprojection::boolean,$)
    global PolySets:
    local Fset,n,i, numlays,R,Rec,CAD1,CAD1Full,CAD1Rec,RecCAD,CurrCADFull,CurrCADRec,cell,cellstack,CADLFull,CADLRec,CADlayer,switch:
    description "LCADRecursiveLR: Computes a recursive layered CAD in listwithrep format. ",
                "Input: A set of polynomials, list of ordered variables, list of cells a previous recursive call terminated at (possibly empty), a layered CAD in listwithrep format from a previous recursive call (possibly empty). Also, an optional parameter stating whether to reset the projection polynomials.",
                "Output: A layered CAD in listwithrep format of one more layer than inputted (or 1 layer if an empty CAD is given), an unevaluated recursive call to produce the following layer (using the 'value' call).":
    
    if kernelopts('opaquemodules')=true then 
        switch:=0: 
        kernelopts('opaquemodules'=false): 
    fi:
    Fset:=convert(F,set):
    #Check if projection polynomials have already been defined
    if type(PolySets,list) then
      userinfo(3, 'LCADRecursiveLR', "Projection polynomials already defined.");
      #If the projection polys have been defined then we need to check the resetproj boolean to see if
      # they need recalculated
      #  if convert(F,set) subset `union`(seq(convert(PolySets[i],set),i=1..nops(PolySets))) then
      if not resetprojection then
        userinfo(3, 'LCADRecursiveLR', "Projection polynomials left as defined.");
      else
        userinfo(3, 'LCADRecursiveLR', "Projection polynomials to be reset as resetproj defined to be true.");
        PolySets:=LCAD_GenerateSeparateProj(Fset,vars,meth):
      end if:
    else
    #If none defined, recompute using LCAD_GenerateSeparateProj
      userinfo(3, 'LCADRecursiveLR', "No projection polynomials previously defined. Recalculating.");
      PolySets:=LCAD_GenerateSeparateProj(Fset,vars,meth):
    end if:
    n:=nops(vars):
    #Calculate the number of layers required
    if nops(C)=0 then
      numlays:=1:
    else
    #Compute numlays by looking at the first cell index and counting the
    # number of even indices, and adding 1.
      numlays:=1+add((i+1) mod 2,i in C[1][1]):
    fi:
    userinfo(3, 'LCADRecursiveLR', "calculated that a", numlays, "layered CAD should be produced.");
    #Create polynomial rings for each layer and split the truncated cell
    # list according to the dimension of each cell
    for i from 1 to n do
      R[i]:=RegularChains:-PolynomialRing(vars[-i..-1]):
      Rec[i]:=select(a->nops(a[1])=i,C):
    end do:
    #Treat the one layered case separately as we have to do the initial
    # splitting of R^1 and dimension checking is slightly simpler.
    if numlays=1 then
      #Generate a full CAD of one-dimensional space (i.e. split R^1) using
      # McCallum's projection method
      CAD1:=CADFull(PolySets[1],[vars[-1]],'method'=meth,'output'='listwithrep'):
      #Separate into full dimensional and deficient cells - store these in a
      # temporary variable and transfer over
      CAD1Full:=select(a->LCAD_IsCellRepOneLayeredDim(a),CAD1);
      CAD1Rec:=remove(a->LCAD_IsCellRepOneLayeredDim(a),CAD1);
      RecCAD:=CAD1Rec:
      CurrCADFull:=CAD1Full:
      CurrCADRec:=CAD1Rec:
      userinfo(3, 'LCADRecursiveLR',"produced", nops(CurrCADFull)+nops(CurrCADRec), "cells of", R[1]['variables'], "-space."):
      userinfo(3, 'LCADRecursiveLR',"discarded", nops(CurrCADRec), "deficient cells of", R[1]['variables'], "-space."):
      #We have to deal separately with the case i=2 as the output of
      # CADFull (used in the separate case above
      for i from 2 to min(2,n) do
        CADLFull:=[]:
        CADLRec:=[]:
        for cell in CurrCADFull do
          #Take each cell in the current list of full cells. Build a stack
          # and separate into full and deficient cells.
          cellstack:=CADGenerateStack(cell,PolySets[i],vars[-i..-1],'output'='listwithrep');
          CADLFull:=[op(CADLFull),op(select(a->LCAD_IsCellRepOneLayeredDim(a),cellstack))]:
          CADLRec:=[op(CADLRec),op(remove(a->LCAD_IsCellRepOneLayeredDim(a),cellstack))]:
        end do:
        #Update temporary variables
        CurrCADFull:=CADLFull:
        CurrCADRec:=CADLRec:
         userinfo(3, 'LCADRecursiveLR',"produced", nops(CurrCADFull)+nops(CurrCADRec), "cells of", R[i]['variables'], "-space."):
        userinfo(3, 'LCADRecursiveLR',"discarded", nops(CurrCADRec), "deficient cells of", R[i]['variables'], "-space."):
        #Add all deficient cells to the list of truncated cells to be
        # symbol '%') for the next layer. (Use type in the future)
       RecCAD:=[op(RecCAD),op(CurrCADRec)]:
      end do:
      for i from 3 to n do
        CADLFull:=[]:
        CADLRec:=[]:
        for cell in CurrCADFull do
          #Take each cell in the current list of full cells. Build a stack
          # and separate into full and deficient cells.
          cellstack:=CADGenerateStack(cell,PolySets[i],vars[-i..-1],'output'='listwithrep');
          CADLFull:=[op(CADLFull),op(select(a->LCAD_IsCellRepOneLayeredDim(a),cellstack))]:
          CADLRec:=[op(CADLRec),op(remove(a->LCAD_IsCellRepOneLayeredDim(a),cellstack))]:
        end do:
        #Update temporary variables
        CurrCADFull:=CADLFull:
        CurrCADRec:=CADLRec:
        userinfo(3, 'LCADRecursiveLR',"produced", nops(CurrCADFull)+nops(CurrCADRec), "cells of", R[i]['variables'], "-space."):
        userinfo(3, 'LCADRecursiveLR',"discarded", nops(CurrCADRec), "deficient cells of", R[i]['variables'], "-space."):
        #Add all deficient cells to the list of truncated cells to be
        # returned in the recursive call. (This will be wrpaped up by a
        # type function later in the package)
        RecCAD:=[op(RecCAD),op(CurrCADRec)]:
      end do:
      #Return the layered CAD along with a recursive call (using the value
      # symbol '%') for the next layer.
      return CurrCADFull,%LCADRecursive(F,vars,RecCAD,CurrCADFull,'method'=meth,'output'=outp,'resetproj'=false):
    fi:
    #Now deal with the case where we already have a layered CAD and are
    # building over the top of it.
    #Store the previous layered CAD in the final output CAD variable (to
    # be updated later)
    CADlayer:=LD;
    #We initialise the variables to be empty as they will be populated by
    # the bottom layer truncated cell after the first loop of the
    # following for loop.
    CurrCADFull:=Rec[1]:
    RecCAD:=[]:
    for i from 2 to n do
      CADLFull:=[]:
      CADLRec:=[]:
      for cell in CurrCADFull do
        #For each cell in the current layer generate a stack using the
        # appropriate projection polynomials and sort into cells with high
        # enough dimension and those deficient.
        cellstack:=CADGenerateStack(cell,PolySets[i],vars[-i..-1],'output'='listwithrep');
        CADLFull:=[op(CADLFull),op(select(a->LCAD_IsCellRepNLayeredDim(a,numlays),cellstack))]:
        CADLRec:=[op(CADLRec),op(remove(a->LCAD_IsCellRepNLayeredDim(a,numlays),cellstack))]:
      end do:
      #Update the current CAD level with those cells generated and those
      # truncated cells inputted into the procedure of the correct dimension.
      CurrCADFull:=CADLFull:
      userinfo(3, 'LCADRecursiveLR',"produced", nops(CurrCADFull)+nops(CurrCADRec), "cells of", R[i]['variables'], "-space."):
      CurrCADFull:=[op(CurrCADFull),op(Rec[i])]:
      CurrCADRec:=CADLRec:
      userinfo(3, 'LCADRecursiveLR',"discarded", nops(CurrCADRec), "deficient cells of", R[i]['variables'], "-space."):
      #Add all deficient cells to the list of truncated cells to be
      # returned in the recursive call.
      RecCAD:=[op(RecCAD),op(CurrCADRec)]:
    end do:
    #Combine the layered CAD inputted into the procedure with newly
    # computed cells.
    CADlayer:=[op(CADlayer),op(CurrCADFull)]:
    if switch=0 then 
        kernelopts('opaquemodules'=true): 
    fi:
    #Return the new layered CAD along with a recursive call (using the
    # value symbol '%') for the next layer.
    return CADlayer,%LCADRecursive(F,vars,RecCAD,CADlayer,'method'=meth,'output'=outp,'resetproj'=false):
end proc:

#####################################################################

LTTI_TTICADNLayered:=proc( PHI::list, N,  vars::list, output::symbol, $ )
    local i, t, n, mvar, lvars, phi, f, gs, E, F, A, B, C, BigF, BigC, P, lowercad, cad, cell, FinalCAD, l, tmpCAD:
    description "LTTI_TTICAD: Compute the TTICAD for a list of clauses.",
                "Input: PHI, a list of lists phi[i], each of which represent a clause.  Each phi[i] has two arguments; the first a polynomial (equational constraint) and the second a list of polynomials (from the other constraints).",
                "Also a list of ordered variables and optionally, an ouput choice."
                "Output: A truth-table invariant CAD for PHI or failure if the polynomials are not well oriented.":

    t := nops(PHI):
    n := nops(vars):
    mvar := vars[1]:
    lvars := remove(X->X=mvar,vars):
    phi := table():
    f := table():
    gs := table():
    E := table():
    F := table():
    for i from 1 to t do
        phi[i] := PHI[i]:
        f[i] := expand(phi[i][1]):
        gs[i] := expand(phi[i][2]):
        E[i] := {f[i]}:
        F[i] := PCAD_PrimSet(E[i], mvar):
        F[i]:= PCAD_SFBasis(F[i], mvar):
    od:
    BigF := `union`(seq(F[i], i=1..t)):
    if n=1 then
        return( LTTI_NLayeredUnivariateCase(BigF, N,mvar, output) ):
    fi:
    A := table():
    B := table():
    C := table():
    for i from 1 to t do
        A[i] := {f[i], op(gs[i])}:
        B[i] := PCAD_PrimSet(A[i], mvar):
        C[i] := PCAD_ContSet(A[i], mvar):
    od:
    BigC := `union`(seq(C[i], i=1..t)):
    P := TTI_TTIProjectionOperator(F, B, t, mvar, lvars):
    P := P union  BigC:
    lowercad := LTTI_NLayeredLowerDimCAD(P,N, lvars, output):
    #print(map(x->x[1],lowercad)):
    if lowercad=FAIL then
        print("The reduced projection set is not well oriented.  Hence CADW fails and the TTICAD cannot be produced"):
        return(FAIL):
    fi:
    cad:=table():
    for i from 1 to nops(lowercad) do
        cell := lowercad[i]:
        cad[i] := LTTI_TTINLayeredGenerateStack( cell, t, mvar, lvars, f, gs, B, E, output):
        if cad[i]=FAIL then
            print("An equational constraint is nullified on a cell of positive dimension and the excluded polynomials are not constant.  Hence this implementation of TTICAD fails."):
            return(FAIL):
        fi:
    od:
    FinalCAD := [entries(cad, 'nolist')]:
    tmpCAD := map(X->op(X), FinalCAD):
    FinalCAD:=[]:
    for l in tmpCAD do
      if LCAD_IsCellRepNLayeredDim(l,N) then
        FinalCAD:=[op(FinalCAD),l]:
      fi:
    od:
    if output='list' or output='listwithrep' then
        return(FinalCAD):
    elif output='rootof' then
        return( map(X->X[2], FinalCAD) ):
    elif output='piecewise' then
        return( PCAD_LWRCADtoPWCAD(FinalCAD) ):
    else
        error("unexpected output format.  This is an unexpected error, please report to D.J.Wilson@bath.ac.uk"):
    fi:
end proc:

#####################################################################

LTTI_TTINLayeredGenerateStack := proc( cell::list, t::posint, mvar::symbol, lvars::list(symbol), f::table, gs::table, B::table, E::table, output, $) :: list:
    local j, vars, Stack, SPrc, SPcube, alpha, beta, falpha, L, cellIndex, excl:
    description "LTTI_TTIGenerateStack: Creating the lifitng set and then stack over a cell in the lower dimensional CAD, for an n-layered TTICAD.":

    cellIndex := cell[1]:
    vars:=[mvar, op(lvars)]:
    SPrc := cell[-1][1]:
    SPcube := cell[-1][2]:
    alpha := PCAD_SPtoRootOf( SPrc, SPcube, RegularChains:-PolynomialRing(lvars) ):
    beta := []:
    for j from 1 to nops(cellIndex) do
        if cellIndex[j]::'even' then
            beta:=[op(beta), alpha[j]]:
        fi:
    od:
    L:={}:
    for j from 1 to t do
        falpha := subs(alpha, f[j]):
        falpha := is(falpha=0):
        if falpha=true then
            userinfo(3, {'TTICAD'},"The polynomial ", f[j], " is nullified on the cell ", cell):
            if PCAD_IsCellZeroDim(cellIndex)=true then
                L := L union B[j]:
                L := remove(X -> evalb(is(subs(alpha,X)=0)), L):
                userinfo(3, {'TTICAD'},"The cell is zero-dimensional so we can continue by expanding the lifting set on this cell."):
            else
                excl := TTI_ExclProj_phi(gs[j], vars):
                excl := remove(X->X::'constant', subs(beta, excl)):
                userinfo(3, {'TTICAD'},"The excluded polynomials are ", excl, " on the cell."):
                if excl = {} then
                    L := L union B[j]:
                    L := remove(X -> evalb(is(subs(alpha,X)=0)), L):
                    userinfo(3, {'TTICAD'},"The relevant excluded polynomials are constants, hence the non-equational constraints are delineable on the cell and we can continue by expanding the lifting set on this cell."):
                else
                    userinfo(3, {'TTICAD'},"The cell has positive dimension and the relevant excluded polynomials are non-constants."):
                    return(FAIL):
                fi:
            fi:
        else
            L := L union E[j]:
        fi:
    od:
    L := PCAD_SetFactors(L):
    L := convert(L, 'list'):
    if output=list then
        Stack := PCAD_GenerateStack(cell, L, RegularChains:-PolynomialRing(vars), 'list'):
    else
        Stack := PCAD_GenerateStack(cell, L, RegularChains:-PolynomialRing(vars), 'listwithrep'):
    fi:
    Stack := map(X->X[1],Stack):
end proc:

#####################################################################

LTTI_NLayeredUnivariateCase := proc( F::set(polynom), N, var::symbol, output::symbol, $ )
    local cad,l,tmpcad:
    description "LTTI_UnivariateCase: Dealing with the n-layered tticad in the case where there is only one variable.":

    if output='list' then
        cad := CADFull( convert(F,'list'), var,'method'='McCallum' ):
        #cad:=map(X->X[1], cad):
        tmpcad:=cad:
        cad:=[]:
        for l in tmpcad do
          if LCAD_IsCellRepNLayeredDim(l,N) then
            cad:=[op(cad),l]:
          fi:
        od:
        return(cad):
    else
        cad := CADFull( convert(F,'list'), var,'method'='McCallum' ):
    fi:
    #tmpcad:=map(X->X[1], cad):
    cad:=[]:
    for l in tmpcad do
      if LCAD_IsCellRepNLayeredDim(l,N) then
        cad:=[op(cad),l]:
      fi:
    od:
    if output='listwithrep' then
        return(cad):
    elif output='piecewise' then
        return( PCAD_LWRCADtoPWCAD(cad) ):
    elif output='rootof' then
        return( map(X->X[2], cad) ):
    else error("unexpected output format.  This is an unexpected error, please report to D.J.Wilson@bath.ac.uk"):
    fi:
end proc:

#####################################################################

LTTI_NLayeredLowerDimCAD := proc( P::set(polynom), N,lvars::list(symbol), output::symbol, $ )
    local pset, cad:
    description "LTTI_LowerDimCAD: Constructing the lower dimension CAD using CADW via LayeredCAD.":

    pset := PCAD_McCProjPolys(P, lvars):
    if output='list' then
      cad:=LCAD_ProjCADLiftNLayered( pset,N, lvars, 'McCallum', 'OI', 0, 'list', 'giveFAIL' ):
    else
      cad:=LCAD_ProjCADLiftNLayered( pset,N, lvars, 'McCallum', 'OI', 0, 'listwithrep', 'giveFAIL' ):
    fi:
    cad:
end proc:

#####################################################################

LTTI_Dist := proc(PHI :: {list, set}, vars :: list) :: 'list':
    local FC,i,DL,cind,cdim,n:
    description "LTTI_Dist: Generates a list with the number of cells of each dimension in a full TTICAD using McCallum projection.",
                "Input: A list or set of polynomials and a list of variables.",
                "Output: A list of the number of cells of each dimension.":

    DL:=[0$(nops(vars)+1)]:
    #Compute a full CAD for the polynomials using McCallum's projection operator
    FC:=TTICAD(PHI,vars, 'output'='list'):
    for i from 1 to nops(FC) do
      #Extract the cell index
      cind:=FC[i][1]:
      #Calculate the dimension by summing the parities of the indices
      cdim:=add(n mod 2, n in cind);
      #Add 1 to the appropriate entry (cdim+1) in the list
      DL[cdim+1]:=DL[cdim+1]+1:
    end do:

    return(DL):
end proc:

#####################################################################
 
LTTI_NormDist :=  proc(PHI :: {list, set}, vars :: list) :: 'list':
    local tmpres,tot, n:
    description "LTTI_NormDist: Generates a list with the normalised number of cells of each dimension in a full TTICAD using McCallum projection.",
                "Input: A list or set of polynomials and a list of variables.",
                "Output: A list of the normalised number of cells of each dimension.":
    #Compute the number of cells using CADDistribution
    tmpres:=LTTI_Dist(PHI,vars):
    #Calculate the total cells
    tot:=add(n, n in tmpres):
    #Divide each entry by the total
    return( tmpres/tot );
end proc:

#####################################################################

MCAD_LiftToManifold:=proc(lowCAD,equCon,vars,outp)
local Tab,TList,i,cell,cellstack,j:
description "MCAD_LiftToManifold: Lifting algorithm to build a manifold CAD.":

Tab := table():
i := 1:
for cell in lowCAD do
  cellstack := CADGenerateStack(cell,[equCon],vars,'output'=outp):
  if nops(cellstack)>1 then
    for j from 1 to (nops(cellstack)-1)/2 do
      Tab[i] := cellstack[2*j]:
      i := i+1:
    end do:
  end if:
end do:
TList := convert(Tab,list):
return(TList):
end proc:

#####################################################################

MCAD_ConstructLowDimCAD:=proc(T,vars)
    local Epoly, Apolys, EProj, lowvars, lowCAD:
    description "MCAD_ConstructLowDimCAD: Constructing the lower dim CAD when building a manifold CAD.":

    Epoly := T[1]:
    Apolys := convert(T[2],set):
    EProj := ECCADProjOp({Epoly},Apolys,vars):
    lowvars:=vars[2..nops(vars)];
    lowCAD := CADFull(EProj,lowvars,'method'='McCallum'):
    return lowCAD:
end proc:

#####################################################################

MCAD_ManifoldCAD:=proc(T,vars)
    local Epoly,lowCAD:
    description "MCAD_ManifoldCAD: Building a Manifold CAD.":

    Epoly := T[1]:
    lowCAD := MCAD_ConstructLowDimCAD(T,vars):
    return MCADLiftOverLowCAD(lowCAD,Epoly,vars):
end proc:

#####################################################################
    
MCAD_MTTICAD:=proc(PHI,vars)
    local TProj,TProjVar,lowvars,lowCAD,EpolyProd,i:
    description "MCAD_MTTICAD: Building a Manifold TTICAD.":

    TProj := TTICADProjFactors(PHI,vars):
    TProjVar := TProj minus select(t->(degree(t,vars[1])>0),TProj);
    lowvars:=vars[2..nops(vars)];
    lowCAD :=CADFull(TProjVar,lowvars,'method'='McCallum'):
    EpolyProd := product(PHI[i][1],i=1..nops(PHI));
    if degree(EpolyProd,vars[1])<1 then
        ERROR("currently a Manifold TTICAD can only be constructed for equational constraints with positive degree in %1.", vars[1]);
    fi:
    return MCADLiftOverLowCAD(lowCAD,EpolyProd,vars);
end proc:

#####################################################################

LMCAD_ConstructLowDimLCAD:=proc(T,n,vars)
    local Epoly, Apolys, EProj, lowvars, lowCAD:
    description "LMCAD_ConstructLowDimLCAD: Constructing the lower dim part of a layered manifold CAD.":

    Epoly := T[1]:
    Apolys := convert(T[2],set):
    EProj := ECCADProjOp({Epoly},Apolys,vars):
    lowvars:=vars[2..nops(vars)];
    lowCAD := LCAD(EProj,n,lowvars,'method'='McCallum'):
    return lowCAD:
end proc:

#####################################################################

LMCAD_LMCAD:=proc(T,n,vars)
    local Epoly,lowLCAD:

    Epoly := T[1]:
    lowLCAD := LMCAD_ConstructLowDimLCAD(T,n,vars):
    return MCADLiftOverLowCAD(lowLCAD,Epoly,vars):
end proc:

#####################################################################

LMTTICAD_LMTTICAD:=proc(PHI,n,vars)
    local TProj,TProjVar,lowvars,lowLCAD,EpolyProd, i:
    description "LMTTICAD_LMTTICAD: Building a Layered Manifold TTICAD.":

    TProj := TTICADProjFactors(PHI,vars):
    TProjVar := TProj minus select(t->(degree(t,vars[1])>0),TProj);
    lowvars:=vars[2..nops(vars)];
    lowLCAD :=LCAD(TProjVar,n,lowvars,'method'='McCallum'):
    EpolyProd := product(PHI[i][1],i=1..nops(PHI));
    if degree(EpolyProd,vars[1])<1 then
        ERROR("currently a Layered Manifold TTICAD can only be constructed for equational constraints with positive degree in %1.", vars[1]);
    fi:
    return MCADLiftOverLowCAD(lowLCAD,EpolyProd,vars);
end proc:

#####################################################################

end module:

CADassistant
============

CADassistant is a tool to enable users to use CAD technology as efficiently as possible. 

Users input their problem into CADassistant using a standardised syntax. They then select a variety of options on the use for their problem. Finally, they select the output forms they wish to have, choosing from:
 - Maple (RegularChains package)
 - Maple (ProjectionCAD package)
 - Maple (SynRAC package)
 - QEPCAD
 - Mathematica
 - Reduce (Redlog package)

CADassistant uses a variety of heuristics and choices to help it's decision procedure. You can also entirely force CADassistant to produce a specific kind of CAD by choosing 'manual' mode. 

Choices:
 - Specific projection operators
 - Standard or partialCAD
 - Standard, single equational constraint, mutltiple equational constraints (implicit), TTICAD.
 - Standard or layered CAD
 - Standard or manifold CAD (if equational constraint present)
 - Quantifier elimination

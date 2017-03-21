
This conclusions were come to when evaluating BISON heavy test
ifba_he_production.base_ifba_only_smeared Zrb2:

* MooseVariable::computeElemValues and friends all have many redundant if
  branches inside loops that frequently result in several unnecessary
  condition checks for inside nested loops (frequently order 10 to 100
  iterations). Fix this by creating a fast-path function templated on the
  boolean variables for all the if branches with instantiations for common
  combinations of the boolean vars.

* linking tcmalloc into libmesh and moose has no measurable impact on
  performance.

* 


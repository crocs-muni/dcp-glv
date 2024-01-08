Scripts for dowloading, parsing and unrolling formulas. It creates a json file for each curve form and operation

1. Clone database:
	`git clone <anonymized for review>`
2. Unroll formulas. It creates folder "unrolled" and a file "formulas_list.txt":
	`python3 unroll.py efd-dcp`
3. Create a list of sources for individual formulas into sources.txt:
	`python3 sources.py efd-dcp`

Unrolled formulas are stored in jsons {form}_{operation}.json in the following dictionary:

{
	coordinate_system:{
		formula_name:{"parameters":list[str], "variables":list[str],"formula":list[string,string]}
...
}

Where formula_name ("projective:madd-2007-bl-3") consists of coordinate_system ("projective") and identifier from efd ("madd-2007-bl-3").
"parameters" and "variables" contain a list of curve parameters (e.g. "a","b" for weierstrass) and variables for the formula (e.g. "X1","Y1").
"formula" list represents left and right side of each line in the formula.


Notes:
   - formula_info.txt contains information about transformations of formulas to a simple affine form (also see README of the efd repository).
   - supported forms: edwards, montgom, shortw, twisted edwards
   - supported operations: doubling, ladder, differential addition, addition
   - some formulas from efd are skipped (broken):

   	shortw jacobian-3 add-1998-hnm
   	shortw jacobian-3 dbl-1998-hnm-2
   	shortw jacobian-3 dbl-1998-hnm
   	edwards projective add-20080225-hwcd
   	shortw jacobian dbl-1998-hnm
   	edwards projective madd-20080225-hwcd
   	shortw jacobian add-1998-hnm
   	shortw jacobian-0 add-1998-hnm
   	shortw jacobian-0 dbl-1998-hnm
   	shortw jacobian-3 dbl-2004-hmv
   	edwards projective add-2007-bl-4
   	edwards projective add-20090311-hwcd

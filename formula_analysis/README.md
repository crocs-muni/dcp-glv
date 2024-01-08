# Classification of formulas

Automatic analysis of curve operations formula with respect to the DCP problem.



### Get formulas

See `formulas/README.md`

### Process formulas

See `process_formulas` notebook which implements the following:

1. For each formula in `unrolling/unrolled` files, we transform the defining rational maps into a unified affine form (see `efd` repository for more details). All transformations are implemented in `utils/formula.py`. Correctness of transformations can be checked by `python -m pytest utils/`. In the notebook, the transformed formulas are loaded using `load_formulas(path)` into a dictionary.
2. For each formula, an object `FormulaZero` is instantiated which collects all polynomials corresponding to intermediate values. Polynomials are then classified as **Hard** or **Easy**:
    - **Easy** polynomials are defined over coordinates of only one point.
    - **Hard** polynomials are the rest. These do not exist for doubling formulas.
    - Some polynomials are classified neither as easy or as hard and are internally classified as `Coordinate` polynomials. These invoke zero in coordinates of any output point in the computation (RPA attack). Since the success of RPA attacks does not depend on the formula but rather the existence of such points on the curve/subgroup, we discard these polynomials. 

    In the notebook, `FormulaZero` objects are loaded using `load_formula0s` and are stored in `formula0s`.

3. For each model and an operation, we filter out (manually) polynomials that are not relevant (e.g., the ZVP does not have solution for nontrivial reasons).
4. `FormulaZero` objects are then stored in `results` directory as dictionaries. 
5. Classification of formulas with respect to ZVP is printed in the notebook. Formulas are divided into: 
    - Resistant (no Hard or Easy polynomials)
    - Vulnerable (containts at least one Easy polynomial),
    - Semi-resistant (contains some Hard polynomials and no Easy).


### Get curves

Run `curves/load_curves.py` to download curves from `dissect.crocs.fi.muni.cz`. Curves are sorted and saved to `curves.json`. 

### Classify formulas and curves

Run `classify_formulas_curves.py`. For each curve form (edwards, shortw, twisted, montgom) and operation (ladd, dadd, add, dbl) a dataframe is stored into `results/form_operation_database.pkl` . Each row in this dataframe contains curve, formula, polynomial and a difficulty entry where:

- `curve` is a name of a standard curve in the given form,
- `formula` is a name of the formula that supports corresponding curve,
- `polynomial` is a polynomial from `FormulaZero` object,
- `difficulty` is either Hard, Easy, or Unsolvable:
  - Hard means that the polynomial is Hard (see above).
  - Easy means that the polynomial is Easy (see above).
  - Unsolvable means that either the polynomial is Easy but with no solution on the curve.

We then classify combinations of formulas and curves as:
    - Resistant - all polynomials (if any) are Unsolvable
    - Vulnerable - contains at least one Easy polynomial,
    - Semi-resistant - contains some Hard polynomials and no Easy ones. This category does not exist for doubling formulas.

Results of this can be seen in `curves_and_formulas notebook`. Two functions are prepared: `recommendation_of_curve` and `recommendation_of_formula` that take curve/formula on input and a dataframe. It outputs "Resistant", "Vulnerable", "Semi-resistant" curves/formulas. For `recommendation_of_curve` it outputs:

- "Vulnerable curves" - a list of curves with which the formula is vulnerable
- "Semi-resistant curves" - a list of curves with which the formula is semi-resistant. This category does not exist for doubling formulas.
- "Resistant curves" - a list of curves that make the given formula resistant











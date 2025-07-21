# noble-tools
A collection of programs which completes the enumeration of the noble polyhedra and assists in their visualization. In addition to the standard Python libraries, the enumeration relies on the NumPy and SymPy libraries. This is mainly programmed in Python, but in a few situations Mathematica is used instead.

## Usage
The following programs may be run if you want to test the enumeration yourself:
* To enumerate all noble polyhedra in orbit types with 0 degrees of freedom, run the following file:
  * `enumerate0D.py`
* To enumerate all noble polyhedra in orbit types with 1 degree of freedom, run the following files in order:
  * `enumerate1D-A.py`
  * `enumerate1D-B.wls`
  * `enumerate1D-C.py`
* The enumeration for orbit types with 2 degrees of freedom is not yet complete. At the moment, you may validate that there are no noble polyhedra in nonmaximal critical equivalence classes by running the following files in order:
  * `initialize2D-A.py`
  * `initialize2D-B.wls`
  * `initialize2D-C.py`
  * `enumerate2D-A.py`

All files and other relevant data will be placed in the `3dmodels` folder when exported. Noble polyhedra are exported in the .OFF format, and additional data on the minimal polynomials of the locations of the orbits of these noble polyhedra will be exported into a `summary.txt` file.

## Additional resources
Additionally, there is a `library` folder containing .OFF files for all noble polyhedra, as well as the data for the minimal polynomials of their locations.
In the future there are plans to add files allowing the conversion of 3D models to other file formats and the generation of figures.

# noble-tools
A collection of programs which completes the enumeration of the noble polyhedra and assists in their visualization. In addition to the standard Python libraries, the enumeration relies on the NumPy and SymPy libraries.

## Usage
The following programs may be run if you want to test the enumeration yourself:
* `enumerate0D.py` - This program fully enumerates all noble polyhedra in orbit types with 0 degrees of freedom.
* `enumerate1D.py` - This program fully enumerates all noble polyhedra in orbit types with 1 degree of freedom.
All files and other relevant data will be placed in the `3dmodels` folder when exported. Noble polyhedra are exported in the .OFF format, and additional data on the minimal polynomials of the locations of the orbits of these noble polyhedra will be exported into a `summary.txt` file.

## Additional resources
Additionally, there is a `library` folder containing .OFF files for all noble polyhedra, as well as the data for the minimal polynomials of their locations.

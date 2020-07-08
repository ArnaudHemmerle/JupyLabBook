# Release notes

### 08/07/2020

FrontendFunctions.py: 0.13
- Add possibility to print absorbers of a scan. 

CustomFunctions.py: 0.13
- Add possibility to print absorbers of a scan.

### 06/07/2020

FrontendFunctions.py: 0.12
- Cosmetic changes in widgets.
- Fix missing '\n' in script import.
- Change the way to import scripts.
- Add possibility to add scan numbers in scripts automatically.
- Add new command 'Insert commands' to extract commands from a log file.

CustomFunctions.py: 0.12
- In XRF: the first spectrum in 'First/last spectrum' plot is now the first non-zero spectrum.


### 30/06/2020

FrontendFunctions.py: 0.11

- Add delay on cell creation to avoid creation of several identical cells on fast computers. 

### 23/06/2020

FrontendFunctions.py: 0.10  
CustomFunctions.py: 0.11

- Add the choice to do a PyNexus fast extract for the extraction of XRF and Isotherms. All the other extractions are set to fast=True by default.
- Add an informative error message when elements are not right in XRF. Set default to 4 instead of 0, 1, 2, 3.
- Add the possibility to print wm ('where motors') through a button. Goes as a table in the pdf.


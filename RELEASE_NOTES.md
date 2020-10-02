# Release notes

### 02/10/2020

CustomFunctions.py: 1.0.2

Corrections for GISX:
- Change the sign of gamma in the calculation of the position of the direct beam after rotation of the detector.
- Always print the positions of the detector in GISX


### 25/09/2020

CustomFunctions.py: 1.0.1
FrontendFunctions.py: 1.0.1

Corrections for XRF:
- Peak identification is now working with more than one element
- Peak are ordered in the caption
- The size of the sheet increases dynamically (previous limit was at 10 peaks)


### 29/07/2020

Release of the version 1.0. From now on, releases will be more sparse and priority will be given to minimization of changes between versions.

CustomFunctions.py: 1.0
FrontendFunctions.py: 1.0

- A user manual is now available in the repo.
- Few typos/display corrections.

### 20/07/2020

CustomFunctions.py: 0.21
- Major update in GIXS/GIXD: keep only the approximate qxy, leave all other options for analysis.
- Automatic zoom in GIXS/Plot_pilatus when the x or y range is wrong.

FrontendFunctions.py: 0.21
- Major update in GIXS/GIXD: keep only the approximate qxy, leave all other options for analysis.
- Possibility not to add scan numbers in print_scripts by leaving the box empty.


### 17/07/2020

JupyLabBook.ipynb:
- Add a test in the first cell to avoid erasing all the info when running it.


### 16/07/2020

FrontendFunctions.py: 0.20
- Insert command now in reversed alphabetic order.
- Removed a bug in XRF when the user remove a peak.
- Show the peaks in the widget XRF when user clicks on validate peaks.
- XRF: User can plot a peak or not with y/n. Possibility to identify peaks with multiple scans.

### 13/07/2020

CustomFunctions.py: 0.20
- Quick fix: the -2 were not removed before profile integration in GIXS and Extract_pilatus. They are now set to 0 before integration (but remains at -2 on the 2D images).

### 12/07/2020

CustomFunctions.py: 0.19
- Removed the binning in GIXS (too complicated for the user, left for analysis).
- Save profiles with header when plotting qxy/q in GIXS.
- Added header in profiles in Plot_pilatus.

FrontendFunctions.py: 0.19
- Removed the binning in GIXS (too complicated for the user, left for analysis).
- Add the possibility to extract the logs into a human-readable format with the command Convert logs.


### 10/07/2020

JupyLabBook.ipynb
- Put a javascript command in the first cell to avoid cell collapse (no need to click on a cell to expand it)

CustomFunctions.py: 0.18
- Add possibility to print peak positions in XRF

FrontendFunctions.py: 0.18
- Add possibility to print peak positions in XRF
- Add next action button when no scan is selected (by mistake)

### 09/07/2020

CustomFunctions.py: 0.17
- Save the profiles in Plot pilatus
- Profile in x and y in Plot_pilatus
- Replaced the pixel image by profiles in Extract_GIXS

FrontendFunctions.py: 0.17
- Fix the display bug in the print_command and add the date of each command
- Profile in x and y in Plot_pilatus
- Replaced the pixel image by profiles in Extract_GIXS

### 08/07/2020

FrontendFunctions.py: 0.16
CustomFunctions.py: 0.16
- Add the possibility to put xmin and xmax in the profile in Extract_Pilatus.

### 08/07/2020

FrontendFunctions.py: 0.15
CustomFunctions.py: 0.15
- Add the possibility to force the value of gamma and delta in GIXS.

### 08/07/2020

FrontendFunctions.py: 0.14
CustomFunctions.py: 0.14
- Save all data in the working directory.

### 08/07/2020

FrontendFunctions.py: 0.13
- Add possibility to print absorbers of a scan.

CustomFunctions.py: 0.13
- Add possibility to print absorbers of a scan.
- Print in red in Vineyard.

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


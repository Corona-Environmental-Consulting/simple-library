# Purpose
This package includes useful functions used for graphing and formatting plots in the Corona style.
The acid/base equilibrium functions are also included.  These functions include unit conversions, pH and alkalinity calculations with chemical addtion,
and calculations of common scaling and corrosivity parameters.

# Initial Installation
Before this package can install correctly, package dependencies must be installed.

## RTools and devtools
On Windows, RTools will have to be installed to build any package from the source.  Follow the instructions on this website.  
Notice that you have to install the software, then add it to your path.

https://cran.r-project.org/bin/windows/Rtools/

The devtools package should be installed to allow for package install from github.  
This package is also useful for making edits to the Corona package.  Run the following line in the RStudio console:

`install.packages("devtools")`

## phreeqc and tidyphreeqc
These packages are used in the CCPP calculation, but neither of them are available on CRAN.  The coronaenv package install won't work without them.

An old version of the phreeqc package must be installed.  Use the following line in the RStudio console:

`install.packages("https://cran.r-project.org/src/contrib/Archive/phreeqc/phreeqc_3.6.3.tar.gz", repos=NULL, type="source")`

Tidyphreeqc must be installed directly from github.  
Ensure that you have the devtools package installed, then run the following line in your RStudio console:

`devtools::install_github("https://github.com/paleolimbot/tidyphreeqc")`

If you want to use the CCPP functionality, you will also need to install the PHREEQC software from USGS.

https://www.usgs.gov/software/phreeqc-version-3

# Package Install
To install the coronaenv package, you must clone this gitlab repository into your Documents folder.  Then, you can install the pacakge by running
the following code in your RStudio console:

`install.packages("~/corona-custom-library", repos = NULL, type = "source")`

# Updates
As the package is revised, you can install an updated version at any time by pulling the code and rerunning the install function.  
Dependencies don't need to be reinstalled unless you are starting with a fresh R install.  

If you want to make your own updates to the package, this repo can be branched and revised with typical git procedures.  
Note that whichever branch you have checked out will be the one that installs when you run the code above.




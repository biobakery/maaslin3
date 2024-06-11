

##### [Replace the contents of this README.MD file with the appropriate "Users manual" needed for the tool]  
  
   
   
# Bioconductor Package Template

## Introduction

> This page gives details concerning guiding principles and formatting required for Bioconductor packages. Also see [Official Bioconductor Package Guidelines](https://bioconductor.org/developers/package-guidelines/) for more information. 
## Requirements

- R-development version can be downloaded using the links below: 
    - [Source](https://stat.ethz.ch/R/daily/)
    - [macOS](https://mac.r-project.org/)
    - [Windows](https://cran.r-project.org/bin/windows/base/rdevel.html)
    
- #### R Studio (recommended) 
	- Direct to the link [HERE](https://rstudio.com/products/rstudio/download/) to get the suitable version of R studio based on the your Operating system. 
	- After the installation is complete, configure a couple of things so that your code will be formatted the way that we prefer it for Bioconductor. 
		-  Set the column width marker to 80 columns. You can find this in the ‘Tools’ menu if you select the ‘Global Options…’ and then look at the ‘Code Editing’ panel. Then make sure that you click the ‘show margin’ option and that it is set to 80 columns. This will help you see if your lines are too long.
    
		-   Second set up the tab to be 4 spaces. You can do this right above where you set the column width marker for 80.   
     
## Installation

> Getting started with the Bioconductor Template: 
- Click on "**Use this Template**" button and fill out the new repository informations (name/description/access). 
- Finally click on "**Create repository from template**" and you now have a new repository based on the Bioconductor Template following the standard layouts. 

Alternatively, 
- Direct to [Github Create New Repository](https://github.com/organizations/biobakery/repositories/new) and select the "**biobakery/Biocondutor-package-template**". Fill out other information (Repository name/description/access)
- Finally click on "**Create repository**" and you now have a new repository based on the Bioconductor Template. 

**NOTE**: Make your repository "**Private**" unless it is ready to be released.

## Developing New Bioconductor Packages 
### Setting up local development environment
- Clone your recently created repository in your local development environment either using:
    ``` 
        git clone https://github.com/biobakery/<Name of your repository>
    ```
    or using the "**Clone or Download**" button. 

### Installing the dependencies
- Get the latest or compatible version of the "**BiocManager**" and "**devtools**". Run the following installation code in the R console using R studio.  
  Current Development Version of Bioconductor: 3.11
    ```
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
          install.packages("devtools"); library("devtools");
        BiocManager::install(version = "3.10")
        BiocManager::valid()              # checks for out of date packages
    ```
    Alternatively, 
    Run the **/R/template.R** file which will install the latest version of "**BiocManager**" and "**devtools**" for you.
- If the project has additional dependencies needed for the development, it is reccommended to be add them to "**Imports**" of the **DESCRIPTION** file. 

## Working with existing packages 
- Detailed instruction to use/clone the Hutlab's bioconductor packages are provided in landing page:
    https://huttenhower.sph.harvard.edu/[toolName]
- Example: 
    ```
    Getting up and running with MaAsLin2 using the bioconductor Template:
    - RUN the template.R file which will install BiocManager and devtools for you. 
    - RUN BiocManager::install("Maaslin2")
    - RUN browseVignettes("Maaslin2") - to view the manual. 
    This will give you the released version of Maaslin2 and you will now be able to run Maaslin2. 
Additional details on the folder stucture and details available [HERE](http://bioconductor.org/developers/how-to/buildingPackagesForBioc/#basic-package-anatomy). 

## Test workflow 
 Two of the main frameworks for testing are [RUnit](http://smarden.org/runit/) and [testthat](https://testthat.r-lib.org/). Additional examples and explanations are provided [here](http://bioconductor.org/developers/how-to/unitTesting-guidelines/). 
 
 For this template, we are using testthat: 

- Once you’re set up the workflow is simple. Modify your code or tests.
- Test your package with ```devtools::test()```
- (Optionally) Test you package using  **Ctrl/Cmd + Shift + T**.  

Repeat until all tests pass.
 - Click "Run Tests" to see a sample demo. 
For more information on testing, see [Offcial R testthat documentation](https://cran.r-project.org/web/packages/testthat/testthat.pdf) or [bioconductor guidelines for Unit test](https://bioconductor.org/developers/package-guidelines/#unittest). 


# Construction and analysis of RNA-seq co-expression networks and development of an interactive interface.
Master's Degree in Big Data and Data Science: Data Science and Data Engineering
Universidad Autónoma de Madrid + IBM

Master Thesis in Google Docs (spanish):

https://docs.google.com/document/d/e/2PACX-1vSqTswCVJBQRs9E3P8s4huU4HrfDk-LjRlC2UxSKMhoGItPgDTYCufXbBDl-SXRkPriQmC44NTN4A-S/pub

Presentation:

https://docs.google.com/presentation/d/1gDQPaApmWmV7bhzEE9IjVv7OYwFDu9p52fLFUj-TfhI/

Feel free to translate it.
...and let me know!

## Using the Shiny app
Currently this Shiny app can be used standalone or in a docker container.
### Standalone
1. Download ``app.R``
2. Download the Rda matrix data for 12 tumors in the link below. You only need to place the Rda files in the same directory with ``app.R`` https://drive.google.com/open?id=1NSIh4fqvXm5bWx7GBqIOGyC4s8dYhiMk
3. Install R
4. You need to install shiny, ape, bioconductor and WGCNA packages. 
5. Packages dplyr and RTCGA are needed if you plan to use the provided data retrieval R script. They will be required for the next version.
6. RStudio is recommended.
7. Run the script either from the command line (``Rscript``), inside R or from Rstudio.
### Docker image
1. Clone repository
2. Download the Rda matrix data for 12 tumors in the link below. You only need to place the Rda files in the same directory with ``app.R`` https://drive.google.com/open?id=1NSIh4fqvXm5bWx7GBqIOGyC4s8dYhiMk
3. ``sudo docker build -t shiny-app .``
4. ``sudo docker run --rm -p 3838:3838 --name shiny-app-deploy shiny-app``
5. Open http://localhost:3838 in your favorite browser.
# Future work
* Docker image should be smaller. I know. It is possible. However, it was not an easy task to prepare a minimum linux with all the needed packages, having to do a lot of manual work. In any case, the tumor data matrices take up too much disk: as a minimum you need more than one gigabyte like it or not. So for the time being one or two GB will do the same pains.
* ETL Data Integration Process to load all TCGA data and keep it up to date. I wanted to dockerize the app first. I haven't made my mind up yet whether the ETL layer will be part of the shiny app or will be an independent process (bet on the last).
* Interface improvements. It is currently a well thought-out interface but only to get out of the way.
  * AND move to Shiny Dashboard! (there may be some changes before going to Shiny Dashboard).
* New functionalities:
  * Add association analysis (see WGCNA documentation and tutorials).
  * Add network and subnetwork (modules) visualization (VisANT and the like).
  * Add differential network analysis (compare two modules). (1)
  * Add network similarity and homology modeling (add evolutionary perspective to co-expression network analysis). This may be an independent but integrated project. (2)

(1) (2) I don't know of any package that solves either : I'm all ears! Some insights about why this adds are at least interesting can be found in my Master Thesis.

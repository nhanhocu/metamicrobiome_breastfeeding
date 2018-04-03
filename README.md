This repository contains the R code and data to reproduce all the analyses and results for the project "Effects of exclusive breastfeeding on infant gut microbiota: a meta-analysis across studies and populations". 

It contains:    
 
1. "metamicrobiome_bf_report.R": the R code that generates all results. In Rstudio, just need to click on the button "compile report" or press Ctrl+Shift+K to generate the html file containing all results (after changing "dir" to your directory).   

2. "metamicrobiome_bf_report.html": the compiled report with all results and R code (generated by "metamicrobiome_bf_report.R").   
 
3. "metamicrobiome_bf_analysis.R": the R code for all the analyses to generate intermediate outputs that are used by "metamicrobiome_bf_report.R" to generate the results in the html file.   
 
4. "data" folder: all data and intermediate outputs generated by "metamicrobiome_bf_analysis.R". These data/outputs are used by "metamicrobiome_bf_report.R" to generate all results.   
 
5. "miscfun.microbiome.R": R functions written specifically for this projects. These R functions are used by "metamicrobiome_bf_analysis.R"  to generate intermediate outputs in 'data' folder and by "metamicrobiome_bf_report.R" to generate all results. A more generally applicable version of these functions are available at: https://github.com/nhanhocu/metamicrobiomeR.   
 

In brief, to reproduce all results in the paper, just click on the the button "compile report" after opening the file "metamicrobiome_bf_report.R" and change "dir" to your directory in Rstudio. 


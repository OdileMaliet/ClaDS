This folder contains the results of our analysis from Maliet, Hartig and Morlon (2019) 
for the 42 MCC bird ClaDS from Jetz (2011) that have more than 50 tips. 

The file "birds_MCC_results.zip" contains a Rdata file for each clade, named "clade_name.Rdata". 
Each of those contains an object "tree" of class "phylo", which is the phylogeny used in the inference, 
and a vector "MAPS", with the inferred parameters of ClaDS2. MAPS[1] is the sigma parameter, MAPS[2] the 
alpha parameter, MAPS[3] the turnover rate epsilon, MAPS[4] the initial speciation rate lambda_0, and 
MAPS[-(1:4)] the branch specific speciation rates, in the same order as tree$edge

The files "clade_name.pdf" show the phylogenies colored with the inferred branch-specific speciation rates
for each of the 42 clade with tip labels. These plots can be obtained in R after loading the corresponding
Rdata file with the following code:

library(RPANDA)
plot_ClaDS_phylo(tree, MAPS[-(1:4)],lwd = 5, show.tip.label = T)

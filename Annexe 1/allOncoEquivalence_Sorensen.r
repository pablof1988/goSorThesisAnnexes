library(goSorensen)
source("..\\adjSignifPvals.R")

# Obtaining ENTREZ identifiers for the gene universe of humans. This vector is 
# directly obtained from the genomic annotation package available in Bioconductor.
# To apply the same process to other species, use the corresponding genomic 
# annotation package. 
library(org.Hs.eg.db)
humanEntrezIDs <- keys(x = org.Hs.eg.db, keytype = "ENTREZID")

# package goSorensen includes object "allOncoGeneLists"
data(allOncoGeneLists)
?allOncoGeneLists
# allOncoGeneLists
sapply(allOncoGeneLists, length)

# ===========================================================================================
# THE FOLLOWING CODE SERVES TO ILLUSTRATIVE PURPOSES, MAINLY. JUMP DIRECTLY TO LINE 140 FOR A
# FASTER AND COMPLETE EQUIVALENCE ANALYSIS OF THESE GENE LISTS
# ===========================================================================================

# The generic function "equivTestSorensen" in package "goSorensen" implements the equivalence test based
# on the Sorensen-Dice distance.
# There are many methods of this function, for different classes of objects passed as arguments.
# Providing two gene lists (essentially, to "character" vectors of gene identifiers) as arguments, it returns
# an object of class "equivSDhtest" inheriting from "htest".
# Providing an object of class "list" of length k with each of its elements representing a gene list (i.e., a "list"
# of k "character" vectors), all possible k*(k - 1)/2 pairwise tests are performed and an object of class "equivSDhtestList"
# (in fact, a "list") of "htest" objects is generated.
# There are also methods for data already summarized in form of 2x2 contingency tables of joint enrichment.

# Examples:
# Equivalence test between gene lists 'waldman' and 'atlas', in dataset 'cancerGeneLists', 
# at level 4 of the BP ontology:
##################################################################
# One Equivalence Test from two Lists Using the Normal Distribution
##################################################################
waldman_atlas.BP.4 <- equivTestSorensen(allOncoGeneLists[["waldman"]], 
                                        allOncoGeneLists[["atlas"]],
                                        geneUniverse = humanEntrezIDs, 
                                        orgPackg = "org.Hs.eg.db",
                                        onto = "BP", GOLevel = 4, 
                                        listNames = c("waldman", "atlas"))
########################
# The Accesor Functions.
########################
# The Sorensen Dissimilarity:
getDissimilarity(waldman_atlas.BP.4)
# The Standard error:
getSE(waldman_atlas.BP.4)
# Equivalence test p-value:
getPvalue(waldman_atlas.BP.4)
# The upper bound for the confidence interval:
getUpper(waldman_atlas.BP.4)
# The enrichment contingency table:
getTable(waldman_atlas.BP.4)

########################
# Bootstrap approach:
########################
boot.waldman_atlas.BP.4 <- equivTestSorensen(allOncoGeneLists[["waldman"]], 
                                             allOncoGeneLists[["atlas"]],
                                             boot = TRUE,
                                             geneUniverse = humanEntrezIDs, 
                                             orgPackg = "org.Hs.eg.db",
                                             onto = "BP", GOLevel = 4, 
                                             listNames = c("waldman", "atlas"))
getTable(boot.waldman_atlas.BP.4)
getPvalue(boot.waldman_atlas.BP.4)
getSE(boot.waldman_atlas.BP.4)

########################
# The upgrade function:
########################
# Building the GO terms mutual enrichment contingency table is the slowest step.
# A faster way to obtain boot.waldman_atlas.BP.4 is to upgrade waldman_atlas.BP.4 to the
# bootstrap approach (but without unnecessarily building again the contingency table):
boot.waldman_atlas.BP.4 <- upgrade(waldman_atlas.BP.4, boot = TRUE)
getPvalue(boot.waldman_atlas.BP.4)
# One can also upgrade another parameters:
# An stricter equivalence limit d0 = 0.2857:
waldman_atlas.BP.4_strict <- upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25)) 
waldman_atlas.BP.4_strict
class(waldman_atlas.BP.4_strict)
# The confidence level
upgrade(waldman_atlas.BP.4, d0 = 1/(1 + 2*1.25), conf.level = 0.99)

############################################################
# One Equivalence Test From the Enriched Contingency Table:
############################################################
# First build the contingency table:
ctab.waldman_atlas.BP.4 <- buildEnrichTable(allOncoGeneLists[["waldman"]],
                                            allOncoGeneLists[["atlas"]],
                                            geneUniverse = humanEntrezIDs,
                                            orgPackg = "org.Hs.eg.db",  
                                            onto = "BP", GOLevel = 4,
                                            listNames = c("waldman", "atlas"))
ctab.waldman_atlas.BP.4

# Now compute the equivalence test:
equivTestSorensen(ctab.waldman_atlas.BP.4)

# Package goSorensen was designed under the object-oriented programming paradigm. All operations,
# like 'equivTestSorensen', have method functions to perform the same operation from adequate classes
# of objects, like performing the equivalence test from two gene lists, from scratch, or to perform
# the equivalence test from the enrichment contingency table associated to them.

################################################################
# All Pairwise Equivalence Tests - Comparing More than two Lists.
################################################################
# All pairwise equivalence tests at level 4 of the BP ontology (quite time consuming, you can jump
# this sentence and use the dataset 'BP.4' [data("BP.4")], which is included in the package 'goSorensen')
# BP.4 <- equivTestSorensen(allOncoGeneLists,
#                          geneUniverse = humanEntrezIDs, 
#                          orgPackg = "org.Hs.eg.db",
#                          onto = "BP", GOLevel = 4)
data(BP.4)

# All p-values in vector form:
getPvalue(BP.4)
# The symmetric matrix of all p-values:
getPvalue(BP.4, simplify = FALSE)

# All Sorensen-Dice dissimiliraties:
getDissimilarity(BP.4)
getDissimilarity(BP.4, simplify = FALSE)

# 95% upper limits of all one-sided confidence intervals for the 
# Sorensen-Dice dissimilarity:
getUpper(BP.4)
getUpper(BP.4, simplify = FALSE)

getSE(BP.4)
getSE(BP.4, simplify = FALSE)

getTable(BP.4)

BP.4_strict <- upgrade(BP.4, d0 = 1/(1 + 2*1.25)) # d0 = 0.2857
BP.4_strict
getPvalue(BP.4_strict)

# All pairwise contingency table at level 4 of the BP ontology (quite time consuming, you can jump
# this sentence and use the dataset 'allTabsBP.4' [data("allTabsBP.4")], which is included in the package 'goSorensen')
# allTabsBP.4 <- buildEnrichTable(allOncoGeneLists,
#                                 geneUniverse = humanEntrezIDs,
#                                 orgPackg = "org.Hs.eg.db",
#                                 onto = "BP", GOLevel = 4)
data(allTabsBP.4)
class(allTabsBP.4)

# From this tableList object, one can compute the equivalence test for all the possible
# pair of lists to compare. The result is exactly the same than obtained in BP.4
equivTestSorensen(allTabsBP.4)

############################################################
# One-run Computations for all Ontologies and GO Levels.
############################################################
# (Extremely time consuming, for the same reason as before: building all
# contingency tables. 
# Alternatively, you may use the dataset 'cancerEquivSorensen' directly,
# it is automatically charged with the package 'goSorensen'),
# By default, the tests are iterated over all GO ontologies and for levels 3 to 10:
# cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists, 
#                                             geneUniverse = humanEntrezIDs, 
#                                             orgPackg = "org.Hs.eg.db")
data("cancerEquivSorensen")

# The same but with the bootstrap approach. Even more time consuming and clearly
# unnecessary, see below for a faster approach:
# set.seed(123)
# boot.cancerEquivSorensen <- allEquivTestSorensen(allOncoGeneLists,
#                                                  boot = TRUE,
#                                                  geneUniverse = humanEntrezIDs,
#                                                  orgPackg = "org.Hs.eg.db")

# It takes its time, but it is much more faster:
set.seed(123)
boot.cancerEquivSorensen <- upgrade(cancerEquivSorensen, boot = TRUE)

# # (Also very time consuming.) It is not required to iterate over all ontologies,
# # "allEquivTestSorensen" iterates these procedures over the specified GO ontologies and levels, 
# # e.g., to iterate only over the GO ontologies MF and BP and levels 4, 5, 6 and 7:
specifcsOntosLevels <- allEquivTestSorensen(allOncoGeneLists, 
                                            ontos = c("MF", "BP"), 
                                            GOLevels = 4:7,
                                            geneUniverse = humanEntrezIDs, 
                                            orgPackg = "org.Hs.eg.db")


# ===========================================================================================
# BUILDING ALL CONTINGENCY TABLES OF JOINT ENRICHMENT FOR ALL THREE GO ONTOLOGIES AND FOR
# LEVELS 3 TO 10
# ===========================================================================================

# Quite time consuming. Jump this line and go to "data(allTabs)" to directly have all tables
# under Bioconductor
# allTabs <- allBuildEnrichTable(allOncoGeneLists,
#                                geneUniverse = humanEntrezIDs, 
#                                orgPackg = "org.Hs.eg.db")
data(allTabs)

# Obtainig all the parwaise equivalence test for the three ontologies, for GO levels
# 3 to 10, very fast:
cancerEquivSorensen <- allEquivTestSorensen(allTabs)
# But the previous sentence may be substituted by
data(cancerEquivSorensen)
# to have all these test results under Bioconductor 3.17 (they depend on the previously built tables).

# The bootstrap version of these tests tends to be more exact (less danger of making a type I
# error) and, in any case, tends to be conservative under low enrichment frequencies. So the positive
# results seem more reliable. The number of valid bootstrap replicates over the initially planned (10000),
# is also displayed. Under low table frequencies, some generated bootstrap tables are not adequate for
# Sorensen-Dice computations, but this only induces a conservative tendency in the test.

# Slightly more time consuming:
set.seed(123)
boot.cancerEquivSorensen <- allEquivTestSorensen(allTabs, boot = TRUE)
# Alternatively, load directly these results
data(boot.cancerEquivSorensen)

##############################################
# Correction For Multiple Comparisons.
##############################################
# The above multiple equivalence tests are performed without any adjustment for
# testing multiplicity, which is left to the responsibility of the user. 
# Some suggestions follow on how to perform this task.

# From 21 possible p-values (21 = 7 * (7 - 1) / 2) and after the Holm's adjustment for 
# testing multiplicity, identify those who are <= 0.05:
# Function adjSignifPvals in script "adjSignifPvals.R" returns the adjusted significant p-values
# jointly with the enrichment contingency tables, in order to put these values in an adequate
# context (e.g., those who are not very credible due to low table frequencies):

# To create a text file with all these results:
sink(file = "signifEquivalencies.txt", append = TRUE)

signifPvals_d0_0.4444 <- adjSignifPvals(cancerEquivSorensen)

cat("========================================================\n")
cat("signifPvals_d0_0.4444$BP\n")
cat("========================================================\n")
signifPvals_d0_0.4444$BP
cat("========================================================\n")
cat("signifPvals_d0_0.4444$CC\n")
cat("========================================================\n")
signifPvals_d0_0.4444$CC
cat("========================================================\n")
cat("signifPvals_d0_0.4444$MF\n")
cat("========================================================\n")
signifPvals_d0_0.4444$MF

# For a more restrictive d0 = 0.2857:
signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(cancerEquivSorensen, d0 = 1/(1 + 2*1.25)))

cat("========================================================\n")
cat("signifPvals_d0_0.2857$BP\n")
cat("========================================================\n")
signifPvals_d0_0.2857$BP
cat("========================================================\n")
cat("signifPvals_d0_0.2857$CC\n")
cat("========================================================\n")
signifPvals_d0_0.2857$CC
cat("========================================================\n")
cat("signifPvals_d0_0.2857$MF\n")
cat("========================================================\n")
signifPvals_d0_0.2857$MF


cat("\n\n\n")

cat("********************************************************\n")
cat("               BOOTSTRAP TEST\n")
cat("********************************************************\n")

boot.signifPvals_d0_0.4444 <- adjSignifPvals(boot.cancerEquivSorensen)

cat("========================================================\n")
cat("boot.signifPvals_d0_0.4444$BP\n")
cat("========================================================\n")
boot.signifPvals_d0_0.4444$BP
cat("========================================================\n")
cat("boot.signifPvals_d0_0.4444$CC\n")
cat("========================================================\n")
boot.signifPvals_d0_0.4444$CC
cat("========================================================\n")
cat("boot.signifPvals_d0_0.4444$MF\n")
cat("========================================================\n")
boot.signifPvals_d0_0.4444$MF

# For a more restrictive d0 = 0.2857:
set.seed(123)
boot.signifPvals_d0_0.2857 <- adjSignifPvals(upgrade(boot.cancerEquivSorensen, d0 = 1/(1 + 2*1.25), boot = TRUE))

cat("========================================================\n")
cat("boot.signifPvals_d0_0.2857$BP\n")
cat("========================================================\n")
boot.signifPvals_d0_0.2857$BP
cat("========================================================\n")
cat("boot.signifPvals_d0_0.2857$CC\n")
cat("========================================================\n")
boot.signifPvals_d0_0.2857$CC
cat("========================================================\n")
cat("boot.signifPvals_d0_0.2857$MF\n")
cat("========================================================\n")
boot.signifPvals_d0_0.2857$MF

sink(file = NULL)

# High level of consistency between normal and bootstrap results is a general trend, 
# with greater p-values in the bootstrap case, as expected. For example:
signifPvals_d0_0.2857$BP$`level 5`
boot.signifPvals_d0_0.2857$BP$`level 5`


####################################################
# IRRELEVANCE-THRESHOLD MATRIX OF DISSIMILARITIES:
####################################################
# For an ontology and GO level in specific, the matrix can be computed from the list
# of enrichement contingency tables:
dismatBP4 <- sorenThreshold(allTabs$BP$`level 4`)
dismatBP4

# The same result, but consuming a huge amount of time (not necessary since 
# the list of enrichment contingency tables are available):
# dismatBP4 <- sorenThreshold(allOncoGeneLists, onto = "BP", GOLevel = 4,
#                            geneUniverse = humanEntrezIDs, 
#                            orgPackg = "org.Hs.eg.db")

# Matrix of dissimilarities using the bootstrap distribution:
dismatBP4.boot <- sorenThreshold(allTabs$BP$`level 4`, boot = TRUE)
dismatBP4.boot

# Equivalent but slower:
# dismatBP4.boot <- sorenThreshold(allOncoGeneLists, onto = "BP", GOLevel = 4,
#                                geneUniverse = humanEntrezIDs, 
#                                orgPackg = "org.Hs.eg.db",
#                                boot = TRUE)

# Dendrograms from these equivalence threshold distances:
clust.threshold <- hclustThreshold(dismatBP4, onTheFlyDev = "windows")
clust.threshold.boot <- hclustThreshold(dismatBP4.boot, onTheFlyDev = "windows")

# All pre-built enrichment contingency tables for all three GO ontologies
# and for GO levels 3 to 10:
data(allTabs)
# These data come from code like:
# allTabs <- allBuildEnrichTable(allOncoGeneLists, 
#                                geneUniverse = humanEntrezIDs, orgPackg = "org.Hs.eg.db")

# Computing all equivalence threshold dissimilarities 
# for these GO ontologies and levels:
allDismat <- allSorenThreshold(allTabs, trace = FALSE)
allDismat$BP
allDismat$BP$`level 4`
# This last matrix should be the same as:`
dismatBP4

# Generating all cluster objects (without displaying them) as a class
# "equivClustSorensenList" object:
oncAllclusts <- allHclustThreshold(allDismat, trace = FALSE)

# Generating pdf files with the plots of all these clusters:
goProfiles::equivClust2pdf(pruneClusts(oncAllclusts), "oncAllclusts")
# Function 'pruneClusts' removes all NULL or not displayable as a dendrogram clusters


# "NMA: Network meta-analysis based on multivariate meta-analysis and meta-regression models in R"
#  by Hisashi Noma
#
#  October 23, 2025
###


# 1. Install and load the "NMA" package

pkgCheck <- function(pkg){
	if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
	library(pkg, character.only = TRUE)
}

pkgCheck("NMA")



# 2. Load the illustrative example
#    Original dataset from: https://doi.org/10.1001/archinternmed.2010.427

data(heartfailure)
print(heartfailure)

# An arm-based dataset of a network meta-analysis that compared 7 antihypertensive drug classes (AB, ACE, ARB, BB, CCB, DD, CT) and placebo.
# For details, please see help(heartfailure).



# 3. Setup: Transforming arm-level data to contrast-based summary statistics and making objects for the network meta-analysis

hf2 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="OR",ref="Placebo",data=heartfailure)
hf3 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="RR",ref="Placebo",data=heartfailure)
hf4 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="RD",ref="Placebo",data=heartfailure)

# The effect measure can be specified by "measure", and the reference treatment can be specified by "ref".



# 4. Creating a network plot

netplot(hf3)                                      # Network plot
netplot(hf3,base.lwd=1.5,base.cex=1.5)            # The edge widths and the node sizes can be changed.
netplot(hf3,col="red",bg="red")                   # The color can be changed.
netplot(hf3,text=FALSE)                           # The text can be cancelled.



# 5. Pairwise meta-analysis for all treatment pairs with direct comparisons

pairwise(hf2)            # The pairwise meta-analyses are performed by "rma" and "regtest" functions of "metafor" package.
pairwise(hf3)
pairwise(hf4)



# 6. Network meta-analysis with the consistency model

nma(hf2)                            # measure = "OR" (odds ratio)
nma(hf2, eform=TRUE)                # The outputs can be changed to exponential scale by "eform".

nma(hf3, eform=TRUE)                # measure = "RR" (risk ratio)

nma(hf4)                            # measure = "RD" (risk difference)




# 7. Network meta-regression

nmareg(hf3,z=SBP,treats=3)
nmareg(hf3,z=c(SBP,DBP),treats=3)
nmareg(hf3,z=c(SBP,DBP),treats=4)
nmareg(hf3,z=c(SBP,DBP),treats=c(3,4,6))




# 8. Creating ranking statistics for network meta-analysis (e.g., SUCRA)

nmarank(hf3)
nmarank(hf3, ascending=FALSE)                                 # The order can be inversed by "ascending".





# 9. Creating a league table

nmaleague(hf3)
nmaleague(hf3, eform=TRUE)                                              # The outputs can be changed to exponential scale by "eform".
nmaleague(hf3, eform=TRUE, digits=2)                                    # Number of decimal places can be changed by "digits"
nmaleague(hf3, eform=TRUE, PI=TRUE)                                     # A league table of prediction intervals can be created by setting "PI=TRUE".
nmaleague(hf3, eform=TRUE, out.csv="nmaleague_out.csv")                 # The outputs can be exported to a CSV file by setting "out.csv='filename'".




# 10. Ranked forest plot

nmaforest(hf3)                                                # Ranked forest plot
nmaforest(hf3, col.plot="blue")                               # The color can be changed.
nmaforest(hf3, ascending=FALSE)                               # The order can be inversed by "ascending".

##

of2 <- obj.forest(hf3)                                        # To make a more customized forestplot, users can use "forestplot" package directly. The input objects can be created by "obj.forest" function.

require(forestplot)

# Example (1): Ticks of x-axis can be changed by "xticks".

forestplot(labeltext=of2$labeltext, of2$coef, boxsize = of2$boxsize,
		align=c("l","c","c"), zero = 1,	
		col = fpColors(lines="blue", box="blue"), 
		txt_gp = fpTxtGp(ticks=gpar(cex=1)),
		xticks=c(0.5,1,1.5,2))


# Example (2): Limits of x-axis can be changed by "clip".

forestplot(labeltext=of2$labeltext, of2$coef, boxsize = of2$boxsize,
		align=c("l","c","c"), zero = 1,	
		col = fpColors(lines="blue", box="blue"), 
		txt_gp = fpTxtGp(ticks=gpar(cex=1)),
		xticks=c(0.5,1,1.5,2,3), clip=c(0.5,3))
		
		

# 11. Local inconsistency tests for all poissible closed loops (generalized Bucher's test)

local.ict(hf3)





# 12. Higgins' global inconsistency test based on the design-by-treatment interaction model

global.ict(hf3)





# 13. Q-statistic and its factorization

nmaQ(hf3)





# 14. Noma's side-splitting

sidesplit(hf3)





# 15. Jackson's random inconsistency model

random.icm(hf3)




# 16. Comparison-adjusted funnel plot

nmafunnel(hf3,legends="bottomright")		# Comparison-adjusted funnel plot for placebo-controlled trials





# 17. Contribution weight matrices

nmaweight(hf3)





# 18. Transitivity analysis

transitivity(hf3, SBP)
transitivity(hf3, SBP, yrange=c(100,220))			# Specify the range of y-axis.
transitivity(hf3, DBP)
transitivity(hf3, pubyear)	

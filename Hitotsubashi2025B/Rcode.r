#------------------------------------------------------------------------
#
# 「ネットワークメタアナリシス」
#  一橋大学ソーシャルデータサイエンス学部
#
# 野間久史　（統計数理研究所）
#
# 2025年11月27日
#
#------------------------------------------------------------------------





# 1. Install the package from CRAN (require R >= 3.5.0)

install.packages("NMA")



# 2. Load the "NMA" package

library("NMA")            # load the "NMA" package



# 3. Load the illustrative example
#    Original dataset from: https://doi.org/10.1001/archinternmed.2010.427

data(heartfailure)
print(heartfailure)

# An arm-based dataset of a network meta-analysis that compared 7 antihypertensive drug classes (AB, ACE, ARB, BB, CCB, DD, CT) and placebo.
# For details, please see help(heartfailure).



# 4. Setup: Transforming arm-level data to contrast-based summary statistics and making objects for the network meta-analysis

hf2 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="OR",ref="Placebo",data=heartfailure)
hf3 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="RR",ref="Placebo",data=heartfailure)
hf4 <- setup(study=study,trt=trt,d=d,n=n,z=c(SBP,DBP,pubyear),measure="RD",ref="Placebo",data=heartfailure)

# The effect measure can be specified by "measure", and the reference treatment can be specified by "ref".



# 5. Creating a network plot

netplot(hf2)                                      # Network plot
netplot(hf2,base.lwd=1.5,base.cex=1.5)            # The edge widths and the node sizes can be changed.
netplot(hf2,col="red",bg="red")                   # The color can be changed.
netplot(hf2,text=FALSE)                           # The text can be cancelled.



# 6. Pairwise meta-analysis for all treatment pairs with direct comparisons

pairwise(hf2,method="REML",test="z")            # The pairwise meta-analyses are performed by "rma" and "regtest" functions of "metafor" package.
pairwise(hf3,method="REML",test="z")
pairwise(hf4,method="REML",test="z")



# 7. Network meta-analysis with the consistency model

nma(hf2)                            # measure = "OR" (odds ratio)
nma(hf2, eform=TRUE)                # The outputs can be changed to exponential scale by "eform".

nma(hf3, eform=TRUE)                # measure = "RR" (risk ratio)

nma(hf4)                            # measure = "RD" (risk difference)




# 8. Network meta-regression

nmareg(hf2,z=SBP,treats=3)
nmareg(hf2,z=c(SBP,DBP),treats=3)
nmareg(hf2,z=c(SBP,DBP),treats=4)
nmareg(hf2,z=c(SBP,DBP),treats=c(3,4,6))




# 9. Creating ranking statistics for network meta-analysis (e.g., SUCRA)

nmarank(hf2)
nmarank(hf2, ascending=FALSE)                                 # The order can be inversed by "ascending".





# 10. Creating a league table

nmaleague(hf2)
nmaleague(hf2, eform=TRUE)                                              # The outputs can be changed to exponential scale by "eform".
nmaleague(hf2, eform=TRUE, digits=2)                                    # Number of decimal places can be changed by "digits"
nmaleague(hf2, eform=TRUE, PI=TRUE)                                     # A league table of prediction intervals can be created by setting "PI=TRUE".
nmaleague(hf2, eform=TRUE, out.csv="nmaleague_out.csv")                 # The outputs can be exported to a CSV file by setting "out.csv='filename'".




# 11. Ranked forest plot

nmaforest(hf2)                                                # Ranked forest plot
nmaforest(hf2, col.plot="blue")                               # The color can be changed.
nmaforest(hf2, ascending=FALSE)                               # The order can be inversed by "ascending".

##

of2 <- obj.forest(hf2)                                        # To make a more customized forestplot, users can use "forestplot" package directly. The input objects can be created by "obj.forest" function.

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
		
		

# 12. Local inconsistency tests for all poissible closed loops (generalized Bucher's test)

local.ict(hf2)





# 13. Higgins' global inconsistency test based on the design-by-treatment interaction model

global.ict(hf2)





# 14. Q-statistic and its factorization

nmaQ(hf2)





# 15. Noma's side-splitting

sidesplit(hf2)





# 16. Jackson's random inconsistency model

random.icm(hf2)




# 17. Comparison-adjusted funnel plot

nmafunnel(hf2,legends="bottomright")		# Comparison-adjusted funnel plot for placebo-controlled trials





# 18. Contribution weight matrices

nmaweight(hf2)





# 19. Transitivity analysis

transitivity(hf2, SBP)
transitivity(hf2, SBP, yrange=c(100,220))			# Specify the range of y-axis.
transitivity(hf2, DBP)
transitivity(hf2, pubyear)	




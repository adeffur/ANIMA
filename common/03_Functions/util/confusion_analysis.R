################################################################
## Analysis of the quality of predictions, on the basis of a table of
## confusion.

## Reference for most metrics:
## http://149.170.199.144/multivar/accuracy/ConfusionMatrix.html

## How I load the command on my computer. Change the directory
## according to your local settings.
## setwd('~/motif_discovery_competition_2003/evaluation/')
## source('confusion_analysis.R')

## Description of the parameters
param.desc <- c(
		## Raw counts
		TP="True positive",
		FP="False positive",
		FN="False negative",
		TN="True negative",

		## Marginal sums
		KP="Known Positive",
		KN="Known Negative",
		PP="Predicted Positive",
		PN="Predicted Negative",

		N="Total",
		
		## Weights (costs)
		TPw="True positive weight",
		FPw="False positive weight",
		FNw="False negative weight",
		TNw="True negative weight",

		## Pseudo-weights
		TPp="True positive pseudo-weight",
		FPp="False positive pseudo-weight",
		FNp="False negative pseudo-weight",
		TNp="True negative pseudo-weight",

		## Corrected weights
		TPc="True positive corrected",
		FPc="False positive corrected",
		FNc="False negative corrected",
		TNc="True negative corrected",
		
		## Corrected maringal sums
		KPc="Known Positive corrected",
		KNc="Known Negative corrected",
		PPc="Predicted Positive corrected",
		PNc="Predicted Negative corrected",
		Nc="Total corrected",
		
		## Evaluation parameters
		Prev="Prevalence",
		ODP="Overall Diagnostic Power ",
		CCR="Correct Classification Rate" ,
		Sn="Sensitivity",
		Sp="Specificity",
		FPR="False Positive Rate",
		FNR="False Negative Rate",
		PPV="Positive Predictive Value",
		NPV="Negative Predictive Value",
		Mis="Misclassification Rate",
		Odds="Odds-ratio",
		Kappa="Kappa",
		NMI="NMI n(s)",
		ACP="Average Conditional Probability",
		MCC="Matthews correlation coefficient",
		Acc1="Accuracy (Sn+Sp)/2",
		Acc2="Accuracy (Sn+PPV)/2",
		Hit.noTN= "Hi rate, without counting the TN (to avoid the effect of their large number)"

#		comment = "Comment"
	      )


param.formula <- c(
		   TP= "TP",
		   FP= "FP",
		   FN= "FN",
		   TN= "TN",

		   KP= "TP+FN",
		   KN= "TN+FP",
		   PP= "TP+FP",
		   PN= "FN+FN",

		   N= "TP + FP + FN + TN",

		   TPw= "TPw",
		   FPw= "FPw",
		   FNw= "FNw",
		   TNw= "TNw",

		   TPp= "TPp",
		   FPp= "FPp",
		   TNp= "FNp",
		   FNp= "TNp",

		   TPc= "weights[1]*(TP+pseudo.weights[1])",
		   FPc= "weights[2]*(FP+pseudo.weights[2])",
		   FNc= "weights[3]*(FN+pseudo.weights[3])",
		   TNc= "weights[4]*(TN+pseudo.weights[4])",

		   KPc= "TPc+FNc",
		   KNc= "TNc+FPc",
		   PPc= "TPc+FPc",
		   PNc= "FNc+FNc",

		   Nc= "Nc <- TPc + FPc + FNc + TNc",

		   Prev= "(TPc + FNc)/Nc",
		   ODP= "(FPc + TNc)/Nc",
		   CCR= "(TPc + TNc)/Nc",
		   Sn= "TPc/(TPc + FNc)",
		   Sp= "TNc/(FPc + TNc)",
		   FPR= "FPc/(FPc + TNc)",
		   FNR="FNc/(TPc + FNc)",
		   PPV= "TPc/(TPc + FPc)",
		   NPV= "TNc/(FNc + TNc)",
		   Mis= "(FPc + FNc)/Nc",
		   Odds= "(TPc + TNc)/(FNc + FPc)",
		   Kappa= "((TPc + TNc) - (((TPc + FNc)*(TPc + FPc) + (FPc + TNc)*(FNc + TNc))/Nc))/(Nc - (((TPc + FNc)*(TPc + FPc) + (FPc + TNc)*(FNc + TNc))/Nc))",
		   NMI= "(1 - -TPc*log(TPc)-FPc*log(FPc)-FNc*log(FNc)-TNc*log(TNc)+(TPc+FPc)*log(TPc+FPc)+(FNc+TNc)*log(FNc+TNc))/(Nc*log(Nc) - ((TPc+FNc)*log(TPc+FNc) + (FPc+TNc)*log(FPc+TNc)))",
		   ACP= "0.25*(Sn+ PPV + Sp + NPV)",
		   MCC= "(TPc*TNc - FPc*FNc) / sqrt[ (TPc+FPc)*(TPc+FNc)*(TNc+FPc)*(TNc+FNc)]",
		   Acc1= "(Sn + Sp/2",
		   Acc2= "(Sn + PPV)/2",
		   Hit.noTN="=TPc/(TPc+FPc+FNc)"
		   )


## ##############################################################
## Calculate statistical parameters for evaluating prediction results
## From a confusion table
evaluate <- function(TP, 
		     FP, 
		     FN, 
		     TN, 
		     weights=c(1,1,1,1),
		     pseudo.weights=c(0,0,0,0),
		     comment=""
		     ) {
    result <- vector()

    ## report parameters
    result["TP"] <- TP
    result["FP"] <- FP
    result["FN"] <- FN
    result["TN"] <- TN

    result["KP"] <- TP+FN
    result["KN"] <- TN+FP
    result["PP"] <- TP+FP
    result["PN"] <- FN+FN

    ## weights
    result["TPw"] <- weights[1]
    result["FPw"] <- weights[2]
    result["FNw"] <- weights[3]
    result["TNw"] <- weights[4]

    ## pseudo-weights
    result["TPp"] <- pseudo.weights[1]
    result["FPp"] <- pseudo.weights[2]
    result["FNp"] <- pseudo.weights[3]
    result["TNp"] <- pseudo.weights[4]

    ## total number of cases
    N <- TP + FP + FN + TN
    result["N"] <- N

    ## include weights and pseudo-weights
    TPc <- weights[1]*(TP+pseudo.weights[1])
    FPc <- weights[2]*(FP+pseudo.weights[2])
    FNc <- weights[3]*(FN+pseudo.weights[3])
    TNc <- weights[4]*(TN+pseudo.weights[4])

    result["TPc"] <- TPc
    result["FPc"] <- FPc
    result["FNc"] <- FNc
    result["TNc"] <- TNc

    result["KPc"] <- TPc+FNc
    result["KNc"] <- TNc+FPc
    result["PPc"] <- TPc+FPc
    result["PNc"] <- FNc+FNc

    ## Weighted sum
    Nc <- TPc + FPc + FNc + TNc
    result["Nc"] <- Nc
    
    ## Prevalence 
    result["Prev"] <- (TPc + FNc)/Nc
    
    ## Overall Diagnostic Power 
    result["ODP"] <- (FPc + TNc)/Nc
    
    ## Correct Classification Rate 	
    result["CCR"] <- (TPc + TNc)/Nc
    
    ## Sensitivity 	
    result["Sn"] <- TPc/(TPc + FNc)
    
    ## Specificity 	
    result["Sp"] <- TNc/(FPc + TNc)
    
    ## False Positive Rate 	
    result["FPR"] <- FPc/(FPc + TNc)
    
    ## False Negative Rate 	
    result["FNR"] <- FNc/(TPc + FNc)
    
    ## Positive Predictive Power 	
    result["PPV"] <- TPc/(TPc + FPc)
    
    ## Negative Predictive Power 	
    result["NPV"] <- TNc/(FNc + TNc)
    
    ## Misclassification Rate 	
    result["Mis"] <- (FPc + FNc)/Nc
    
    ## Odds-ratio 	
    result["Odds"] <- (TPc + TNc)/(FNc + FPc)
        
    ## Kappa 	
    result["Kappa"] <- ((TPc + TNc) - (((TPc + FNc)*(TPc + FPc) + (FPc + TNc)*(FNc + TNc))/Nc))/(Nc - (((TPc + FNc)*(TPc + FPc) + (FPc + TNc)*(FNc + TNc))/Nc))
    
    ## NMI n(s) 	
    result["NMI"] <- (1 - -TPc*log(TPc)-FPc*log(FPc)-FNc*log(FNc)-TNc*log(TNc)+(TPc+FPc)*log(TPc+FPc)+(FNc+TNc)*log(FNc+TNc))/(Nc*log(Nc) - ((TPc+FNc)*log(TPc+FNc) + (FPc+TNc)*log(FPc+TNc)))
    
    ## Average Conditional Probability
    result["ACP"] <- 0.25*(result["Sn"] + result["PPV"] + result["Sp"] + result["NPV"])

    ## Matthe's correlation coefficient
    result["MCC"] = (TPc*TNc - FPc*FNc) /((TPc+FPc)*(TPc+FNc)*(TNc+FPc)*(TNc+FNc))^0.5

    ## Accuracy (Sn+Sp)/2
    result["Acc1"] <- (result["Sn"] + result["Sp"])/2

    ## Accuracy (Sn+PPV)/2
    result["Acc2"] <- (result["Sn"] + result["PPV"])/2

    ## Hit rate without counting the TN
    result["Hit.noTN"]=TPc/(TPc+FPc+FNc)

    ## Add comment
#    result["comment"] <- comment

    return (result)
}



################################################################
## Compute the confusion table and hit rate by comparing a vector of
## true classes with a vector of predicted classes.
calc.hit.rate <- function (true.class,pred.class) {
  true.class <- as.vector(true.class)
  pred.class <- as.vector(pred.class)

  if (length(true.class) != length(pred.class)) {
    stop("Error in procedure calc.hit.rate(): true class and pred.class vectors must have the same length.")
  }
  

  result <- list()
  result$confusion.table <- table("true"=true.class,"predicted"=pred.class)
#  validation.hit.rate <- sum(diag(validation.pred.table))/sum(validation.pred.table)
#  print(paste("Validation set hit rate", validation.hit.rate))

  result$n <- length(true.class)
  result$hits <- sum(true.class == pred.class)
  result$errors <- sum(true.class != pred.class)
  result$hit.rate <- result$hits / result$n
  return(result)
}



################################################################
## Quick test
## raw data

L <- 5000
sites <- 10
site.len <- 9
K <- sites*site.len

################################################################
## Test typical behaviours
one.test <- function (weights = c(1,1,1,1),
                      pseudo.weights = c(0,0,0,0)) {

  result <- data.frame(
                       ## ##############################################################
                       ## positive control
                       perfect= evaluate(TP=90, FP=0, FN=0, TN=4910, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="The perfect prediction"),
                       partial= evaluate(TP=50, FP=0, FN=40, TN=4910, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="completely false prediction"),
                       overlap= evaluate(TP=30, FP=60, FN=60, TN=4850, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="completely false prediction"),
                       overpred= evaluate(TP=90, FP=200, FN=0, TN=4710, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="over-prediction encompassing the correct"),
                       nopred= evaluate(TP=0, FP=0, FN=90, TN=4910, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="nopred : no prediction at all"),
                       overpred.overlap= evaluate(TP=30, FP=220, FN=60, TN=4690, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="over-prediction encompassing the correct"),
                       wrong= evaluate(TP=0, FP=90, FN=90, TN=4820, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="completely false prediction"),
                       wrong.overpred= evaluate(TP=0, FP=300, FN=90, TN=4610, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="completely false prediction"),
                       all= evaluate(TP=90, FP=4910, FN=0, TN=0, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="trivial over-predictor: predict everything"),
                       
                       ## ##############################################################
                       ## Negative control: sequences without implants
                       noimp.nopred= evaluate(TP=0, FP=0, FN=0, TN=5000, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="nopred : no prediction at all"),
                       noimp.wrong= evaluate(TP=0, FP=90, FN=0, TN=4910, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="completely false prediction"),
                       noimp.overpred= evaluate(TP=0, FP=300, FN=0, TN=4700, weights=weights, pseudo.weights=pseudo.weights, 						  comment="over-prediction encompassing the correct"),
                       noimp.all= evaluate(TP=0, FP=5000, FN=0, TN=0, weights=weights, pseudo.weights=pseudo.weights, 
                         comment="trivial over-predictor: predict everything")
			 
                       ## ##############################################################
                       ## Description of the parameters
                                        #			 description = param.desc,
                                        #			 formula = param.formula
                       )

  ## transpose the result for readability
  result <- as.data.frame(t(result))
  return (result)
}



do.some.tests <- function () {
  ## export the description and formula
  write.table(data.frame(param.desc,param.formula), file="descriptions.tab", quote=F, sep='\t')

  selected <- c("Sn", "Sp","PPV","NPV","ACP")

  ## Default behaviour: no weight no pseudo-weight
  test0 <- one.test(weights=c(1,1,1,1), pseudo.weights=c(0,0,0,0))
  print("test0"); print(test0[,selected])
  write.table(test0, file="test0.tab", quote=F, sep='\t')

  ## equal pseudo-weights and weights
  test1 <- one.test(weights=c(1,1,1,1), pseudo.weights=c(0.1,0.1,0.1,0.1))
  print("test1"); print(test1[,selected])
  write.table(test1, file="test1.tab", quote=F, sep='\t')

  ## differential weights, no pseudo-weight
  test2 <- one.test(weights=c(10,5,5,1), pseudo.weights=c(0,0,0,0))
  print("test2"); print(test2[,selected])
  write.table(test2, file="test2.tab", quote=F, sep='\t')

  ## differential weight + equal pseudo-weight
  test3 <- one.test(weights=c(10,5,5,1), pseudo.weights=c(0.1,0.1,0.1,0.1))
  print("test3"); print(test3[,selected])
  write.table(test3, file="test3.tab", quote=F, sep='\t')

  ## pseudo weight on TP only
  test4 <- one.test(weights=c(1,1,1,1), pseudo.weights=c(0.1,0,0,0))
  print("test4"); print(test4[,selected])
  write.table(test4, file="test4.tab", quote=F, sep='\t')

  ## differential pseudo weights
  test5 <- one.test(weights=c(1,1,1,1), pseudo.weights=c(0.1,0.1,0.01,0.01))
  print("test5"); print(test5[,selected])
  write.table(test5, file="test5.tab", quote=F, sep='\t')

}

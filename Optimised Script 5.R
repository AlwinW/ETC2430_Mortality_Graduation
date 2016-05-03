#=========================================
#--- Load relevant packages
#-----------------------------------------
# Demography Package for hmd.mx, etc
install.packages("demography")
library("demography")
# RCurl Package
install.packages("RCurl")
library("RCurl")
# fmsb Package for Gompertz Makeham model
install.packages("fmsb")
library("fmsb")
# grDevices Package for exporting png files
install.packages("grDevices")
library("grDevices")

#=========================================
#--- Initial Countries and Variables
#-----------------------------------------
Countries <- matrix(c("Australia", "AUS", "Ireland", "IRL",
                      "Italy", "ITA", "Japan", "JPN",
                      "Russia", "RUS", "Switzerland", "CHE",
                      "U.S.A.", "USA", "Ukraine", "UKR"),
                    nrow = 2, ncol = 8)
# coun = 1
years = 1999:2009
ages = cbind(0:90)

for (coun in 8:ncol(Countries))
{
  #=========================================
  #--- Initial Countries and Variables
  #-----------------------------------------
  # Load mortility data from Mortility.org
  suppressWarnings(rawmortdata <- hmd.mx(Countries[2,coun], "user_email","user_password", Countries[1,coun]))
  
  m = length(years)
  n = length(ages)
  
  # Reduce the data to the required range
  mortdata <- extract.years(rawmortdata,years)
  mortdata <- extract.ages(mortdata,ages,FALSE) # Setting it to True will combine the rates above this range
  m.qx <- mortdata$rate$male
  f.qx <- mortdata$rate$female
  qx <- array(c(m.qx,f.qx) , dim = c(dim(m.qx),2))
  gender <- c("Males", "Females")
  
  # Find a model for each gender
  for (g in 1:2)
  {
    # Initialise the beta and qhat vectors as zero vectors to be filled later
    beta1=(1:m)*0; beta2=(1:m)*0; beta3=(1:m)*0
    qhat.x=array(0,dim=c(n,m))
    
    # Fit the model for each year
    for (i in 1:m)
    {
      model <- fitGM(,qx[1:n,i,g])
      beta1[i] <- model[1]
      beta2[i] <- model[2]
      beta3[i] <- model[3]
      qhat.x[,i] = GompertzMakeham(beta1[i],beta2[i],beta3[i],ages)	
    }
    
    MAPE = sum(abs(qhat.x - qx[,,g])/qx[,,g])/(n*m)
    MaPE = sum((qhat.x - qx[,,g])/qx[,,g])/(n*m)
    
    # Save the data
    # Find years that with reasonable Gompertz Makeham fits
    AC <-  beta1*beta3
#     yearsindex <- 1:m
#     yearsindex <- yearsindex[AC>0]
#     l = yearsindex[1]
#     c = yearsindex[ceiling(length(yearsindex)/2)]
#     u = yearsindex[length(yearsindex)]
    l=1
    c=floor(m/2)
    u=m

    # Set the main title
    plottitle <- paste(sub=Countries[1,coun],"-",gender[g], "aged", head(ages,n=1), "to", tail(ages,n=1), 
                       "in", head(years,n=1), "to", tail(years,n=1), sep=" ")
    # Set the title for GM
    GMtitle <- paste("Gomperz-Makeham model shown for",
                     years[l], ",",years[c], "and", years[u], sep=" ")
    # Set the title for MAPE and MaPE
    subtitle <- paste("MAPE:",MAPE,"and MaPE:",MaPE,sep=" ")
    
    #--- Plot the orginal qx for the data
    png(filename=paste("~/R/",plottitle,".png",sep=""), bg="white")
    
    # Plot the orgininal qx data
    matplot(ages,matrix(log(qx[,,g]),nrow=n,ncol=m), type="l", lty = 1, lwd=1, col=rainbow(m*1.25),
            xlab="Age", ylab="log of Mortality, ln(qx)", main=plottitle, 
            frame=TRUE, axes=TRUE, ylim=c(-10,-1))
    # Create the legend
    legend("topleft",legend=unique(years),
           col=rainbow(m*1.25), ncol=4, pch=19, 
           title="Year", cex=0.75)
    dev.off()
    
    #--- Plot qhat data for three years on the same plot
    png(filename=paste("~/R/",plottitle," - GM",".png",sep=""), bg="white")
    matplot(ages[1:n-1],log(qhat.x[1:n-1,c(l,c,u)]), type="l", lty = 1, lwd=1, col=rainbow(3*1.25),
            xlab="Age", ylab="log of Fitted Mortality, ln(qx)", main=GMtitle, sub=subtitle,
            frame=TRUE, axes=TRUE, ylim=c(-10,-1))
    legend("topleft",legend=years[c(l,c,u)],
           col=rainbow(3*1.25), ncol=3, pch=19, 
           title="Year", cex=0.75)
    dev.off()
    
    # Calculate the Residuals
    residuals = qhat.x - qx[,,g]
    #--- Plotting the residuals with and without the errors
    if (length(AC[AC>0])>0)
    {
    png(filename=paste("~/R/",plottitle," - res smooth only",".png",sep=""), bg="white")
    matplot(ages[5:(n-1)],residuals[5:(n-1),AC>0], pch=1, col=rainbow(length(AC[AC>0])*1.25),
            xlab="Age", ylab="Residual of Fitted Mortality, qhat - qx", main=paste(plottitle,"Residuals",sep=" - "),
            frame=TRUE, axes=TRUE,ylim=c(-0.005,0.005))
    abline(h=0)
    legend("bottomleft",legend=unique(years[AC>0]),
           col=rainbow(length(AC[AC>0])*1.25), ncol=4, pch=19, 
           title="Year", cex=0.65)
    dev.off()
    }
    
    if (length(AC[AC<0])>0)
    {
    png(filename=paste("~/R/",plottitle," - res errors only",".png",sep=""), bg="white")
    matplot(ages[5:(n-1)],residuals[5:(n-1),AC<0], pch=1, col=rainbow(length(AC[AC<0])*1.25),
            xlab="Age", ylab="Residual of Fitted Mortality, qhat - qx", main=paste(plottitle,"Residuals",sep=" - "),
            frame=TRUE, axes=TRUE ,ylim=c(-0.05,0.05))
    abline(h=0)
    legend("bottomleft",legend=unique(years[AC<0]),
           col=rainbow(length(AC[AC<0])*1.25), ncol=9, pch=19, 
           title="Year", cex=0.65)
    dev.off()
    }
    
    # Third Difference
    thirddif  <- qhat.x[4:(n-1),]-3*qhat.x[3:(n-2),]+3*qhat.x[2:(n-3),]-qhat.x[1:(n-4),]
    average  <- (qhat.x[4:(n-1),]+qhat.x[3:(n-2),]+qhat.x[2:(n-3),]+qhat.x[1:(n-4),])/4
    relthirddif <- thirddif/average
    
    # Summarise the data in terms of quantiles
    print(Countries[1,coun])
    print(quantile(relthirddif,probs = c(0.25,0.5,0.75)))
    }
}
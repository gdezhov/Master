library(PowerTOST)

###Sample size for fix sample
#CV = 0.1
theta0 = 0.95
CV = 0.25

Nfix = sampleN.TOST(alpha = 0.05,
				targetpower = 0.8, 
				logscale = TRUE, 
				theta0=theta0, 
				theta1=.8, 
				theta2=1.25, 
				CV=CV, 
				design = "parallel", 
				method="exact", 
				robust=FALSE, 
				print = TRUE, 
				details = FALSE, 
				imax=100)[7]
Nfix
				
#CV = 0.25				
theta0 = 0.95
CV = 0.25



#### Simulation
#gernerate df with 10 000 and 100 N(mu, sigma)
# inputs
##################################
N_sam=100
sim.num=10000
CV=0.25
#Clculations
sd=sqrt(log(1+CV^2))
mu=log(290)
#Log-Normally distributed
X=rlnorm(N_sam,mu,sd)
hist(X)


Y=log(X)
#hist(Y)
shapiro.test(Y)
#########################
#Simulation for Test
#########################
sd=sqrt(log(1+CV^2))
N_sam=1000
sim.num=10000
mu_t = log(290)
mu_r = log(290/0.95)

m <- matrix(0, ncol = N_sam, nrow = sim.num)
dft <- data.frame(m)
#100 000 generated normals
for (i in 1:sim.num) {
set.seed(i)
Xt=rlnorm(N_sam,mu_t,sd)
dft[i,]=Xt
}

path = "C:\\Users\\georgi.dezhov\\Documents\\Dipl\\OUTPUT\\"
out_dft = paste0(path, "dft.csv")
write.csv(dft,out_dft, row.names = FALSE)

#########################
#Simulation for Reference
#########################
m <- matrix(0, ncol = N_sam, nrow = sim.num)
dfr <- data.frame(m)
#100 000 generated normals
for (i in 1:sim.num) {
set.seed(i+10000)
Xr=rlnorm(N_sam,mu_r,sd)
dfr[i,]=Xr
}

out_dfr = paste0(path, "dfr.csv")
write.csv(dfr,out_dfr, row.names = FALSE)

#########################

N_arm = as.numeric(ceiling(Nfix*1.11/2/2))
N_arm
out_t <- matrix(0, ncol = 6, nrow = nrow(df))
out_t <- data.frame(out_t)
colnames(out_t) <- c("GMR", "VAR", "GMR2", "VAR2", "GMR2_SD", "VAR2_SD")
#Geometric Mean of Y (Not LN transformed Data) GMR
out_t$GMR=exp(rowMeans(log(dft[,1:N_arm]))) 
#Var Test
for (i in 1:nrow(df)) {
out_t[i,2]=var(log(as.numeric(dft[i,1:N_arm]))) #Vpool
}

out_r <- matrix(0, ncol = 6, nrow = nrow(df))
out_r <- data.frame(out_r)
colnames(out_r) <- c("GMR", "VAR", "GMR2", "VAR2", "GMR2_SD", "VAR2_SD")
#Geometric Mean of Y (Not LN transformed Data) GMR
out_r$GMR=exp(rowMeans(log(dfr[,1:N_arm]))) 

#Var reference
for (i in 1:nrow(df)) {
out_r[i,2]=var(log(as.numeric(dfr[i,1:N_arm]))) #Vpool
}

out <- matrix(0, ncol = 18, nrow = nrow(df))
out <- data.frame(out)
colnames(out) <- c("GMR", "Vpool", "SDpool", "CVpool", "SE", "POWER", "L", "U", "PASS", "N2",
					"Ntot", "GMRk2", "Vpoolk2", "SEk2", "Lk2", "Uk2", "PASSk2", "OverAll")
out$GMR = out_t[,1]/out_r[,1]

GMR = out$GMR
#####################################
# JPEG device
jpeg("GMR.jpeg", quality = 75)

# Code
hist(GMR)

# Close device
dev.off()
#####################################

out$Vpool = (out_t[,2]+out_r[,2])/2
out$SDpool=sqrt(out[,2])					#SDpool
out$CVpool=sqrt(exp(out[,2])-1)				#CVpool

CV_pooled = out$CVpool
#####################################
# JPEG device
jpeg("CV_pooled.jpeg", quality = 75)

# Code
hist(CV_pooled)

# Close device
dev.off()
#####################################
#Standart Error
out$SE=sqrt(2*out[,2]/N_arm) #SE

###alpha=0.05
alpha=0.05
#Number of observations and degree of freedom
DF=2*N_arm-2
N1=2*N_arm
# t(1-alpha,df)
t_apha_df=qt(1-alpha, DF)

###Power at alpha=0.05
N1 = 2*N_arm
theta2 = 1.25 
theta1 = 0.8
for (j in 1:nrow(out)) {
	out[j,6]=as.numeric(round(pt((log(theta2)-log(out[j,1]))/(out[j,3]*sqrt(2/(N1/2)))-t_apha_df,DF)-
								pt((log(theta1)-log(out[j,1]))/(out[j,3]*sqrt(2/(N1/2)))+t_apha_df,DF),4))
}

for (i in 1:nrow(out)) {

	if (round(as.numeric(out[i,6]),4)>=0.8 & round(out[i,1],2) > 0.8 & round(out[i,1],2) < 1.25) {
	alpha=0.05

	#Lower Bound
	out[i,7] = round(exp(log(out$GMR[i]) - qt(1-alpha, DF)*out[i,5]),2)
					
	#Upper Bound
	out[i,8] = round(exp(log(out$GMR[i]) + qt(1-alpha, DF)*out[i,5]),2)
	
	} else {
	alpha=0.0294

	#Lower Bound
	out[i,7]=round(exp(log(out$GMR[i])- qt(1 - alpha, DF)*out[i,5]),2)

	#Upper Bound
	out[i,8]=round(exp(log(out$GMR[i]) + qt(1 - alpha, DF)*out[i,5]),2)
	
	}
}	

for (i in 1:nrow(out)) {
	
	if (out[i,7] > 0.8 & out[i,8] <1.25 & round(out[j,1],2) > 0.8 & round(out[j,1],2) < 1.25) {
	out[i,9] = as.character("PASS 1")
	} else {
	out[i,9] = as.character("CHECK")
	}
}

for (i in 1:nrow(out)) {
	if (round(out[i,6],4)<0.8 & round(out[i,1],2) > 0.8 
			& round(out[i,1],2) < 1.25 & out[i,9] == ("CHECK")) {
	Nwant=as.numeric(sampleN.TOST(alpha = 0.0294, 
									targetpower = 1-beta, 
									logscale = TRUE, 
									theta0=out[i,1], 
									theta1=theta1, 
									theta2=theta2, 
									CV=out[i,4], 
									design = "parallel", 
									method="exact", 
									robust=FALSE, 
									print = FALSE, 
									details = FALSE, 
									imax=100)[7])
	out[i,10]=Nwant-N1
	out[i,11]=Nwant						#Ntot
	out[i,9] = as.character("Stage 2")
	} else if (out[i,9] == ("CHECK")) {
	out[i,10]=0
	out[i,9] = as.character("FAIL")
	out[i,11]=N1						#Ntot
	} else {
	out[i,10]=0
	out[i,11]=N1
	}
}



### Stage 2
#Geometric Mean of Y (Not LN transformed Data) GMR K2
##more than 1000 needed to be replaced to 1000
TNarm = as.numeric(out[,11]/2)
#length(TNarm[TNarm13000])
TNarm[TNarm>1000]=1000
k=list()
for (i in 1:nrow(out)){
	
	if(out[i,9] == "Stage 2" ){
		k[[i]]=seq(1, TNarm[i], 1)
		
	out_r$GMR2[i] = exp(rowMeans(log(dfr[i,(as.numeric(k[[i]]))]))) 
	out_t$GMR2[i] = exp(rowMeans(log(dft[i,(as.numeric(k[[i]]))]))) 
	out$GMRk2 = out_t[,3]/out_r[,3]

	out_r[i,4]=var(log(as.numeric(dfr[i,(as.numeric(k[[i]]))]))) 		#Vpoolk2 Test 
	out_t[i,4]=var(log(as.numeric(dft[i,(as.numeric(k[[i]]))]))) 		#Vpoolk2 Reference
	out$Vpoolk2 = (out_t[,4]+out_r[,4])/2
	out[i,14]=sqrt(2*out[i,13]/out[i,11]/2) 				#SEk2
	
	alpha=0.0294
	DF=c()
	DF[i]=out[i,11]-2
	out[i,15]=round(exp(log(out$GMRk2[i])-qt(1-alpha, DF[i])*out[i,14]),2)	#Lower Bound
	out[i,16]=round(exp(log(out$GMRk2[i])+qt(1-alpha, DF[i])*out[i,14]),2)	 #Upper Bound
	
	} else {
	out[i,12] = 0	#GMRk2
	out[i,13] = 0	#Vpoolk2
	out[i,14] = 0	#SEk2
	out[i,15] = 0	#Lower Bound
	out[i,16] = 0 	#Upper Bound
	}
}

#"EoStg2"
for (i in 1:nrow(out)) {
	if (out[i,9] == "Stage 2" & (out[i,15] < 0.8 | out[i,16] > 1.25)) {
	out[i,17] = as.character("FAILk2")
	} else if ( out[i,9] == "Stage 2" & (out[i,15] >= 0.8 | out[i,16] <= 1.25)) {
	out[i,17] = as.character("PASSk2")
	} else if (out[i,9] == "FAIL" ) {
	out[i,17] = as.character("FAILk1")
	} else {
	out[i,17] = as.character("PASSk1")
	}
}

#OverAll
for (i in 1:nrow(out)) {
	if (out[i,17] == "FAILk1" |out[i,17] == "FAILk2") {
	out[i,18] = as.character("FAIL")
	} else {
	out[i,18] = as.character("PASS")
	}
}

#Export

output = paste0("C:\\Users\\georgi.dezhov\\Documents\\Dipl\\OUTPUT\\AD__", "_", N_arm, ".csv")
write.csv(out,output, row.names = TRUE)

ntabl = c(min(out[,11]), quantile(out[,11],  probs = 0.05), median(out[,11]),
		mean(out[,11]),quantile(out[,11],  probs = 0.95), max(out[,11]), length(out[,11][out[,11]<=54]))
names(ntabl)=c("Min", "5%", "Median", "Mean", "95%", "Max", "<=54")
ntabl

ntab2 = c(min(out[,10]), quantile(out[,10],  probs = 0.05), median(out[,10]),
		mean(out[,10]),quantile(out[,10],  probs = 0.95), max(out[,10]), length(out[,10][out[,10]<=54]))
names(ntab2)=c("Min", "5%", "Median", "Mean", "95%", "Max", "<=54")
ntab2


ADtabl <- matrix(0, ncol = 7, 2)
ADtabl = data.frame(ADtabl)
ADtabl[2,] = ntabl
ADtabl[1,] = ntab2
rownames(ADtabl) = c("Stage 2", "Total")
colnames(ADtabl) = c("Min", "5%", "Median", "Mean", "95%", "Max", "<=54")
write.csv(ADtabl,"C:\\Users\\georgi.dezhov\\Documents\\Dipl\\OUTPUT\\ADtabl.csv", row.names = TRUE)
#####################################
Total_sample_Adaptive = out[,11]
# JPEG device
jpeg("Total_sample_AD.jpeg", quality = 75)

# Code
hist(Total_sample_Adaptive, breaks = 250, main = "Total sample Adaptive design")

# Close device
dev.off()
#####################################

######################################################################
#Adaptive design
######################################################################
Nfix = as.numeric(Nfix)
gsDesign(k = 2, test.type = 2,sfu= "Pocock", n.fix = Nfix, alpha=0.025,  beta = 0.2)


out_gS <- matrix(0, ncol = 21, nrow = nrow(df))
out_gS <- data.frame(out_gS)
colnames(out_gS) <- c("GMR", "Vpool", "SDpool", "CVpool", "SE", "alpha1", "t_apha_df1", "L", "U", "EoStg1",
					"N2", "Ntot", "GMRk2", "Vpoolk2", "SEk2", "alpha2", "t_apha_df2", "Lk2", "Uk2", "EoStg2", "OverAll")

#Geometric Mean of Y (Not LN transformed Data) GMR
out_gS$GMR = out_t[,1]/out_r[,1]

#Pooled Var, SD and CV
out_gS$Vpool = (out_t[,2]+out_r[,2])/2
out_gS$SDpool=sqrt(out_gS[,2])					#SDpool
out_gS$CVpool=sqrt(exp(out_gS[,2])-1)				#CVpool

#Standart Error
out_gS$SE=sqrt(2*out_gS[,2]/N_arm) #SE

###alpha1=0.0294
DF=2*N_arm-2
alpha1=0.0294
out_gS$alpha1=alpha1

# t(1-alpha1,df)
out_gS$t_apha_df1=qt(1-alpha1, DF)

#Lower Bound
out_gS$L=round(exp(log(out_gS$GMR)-(qt(1-out_gS$alpha1, DF))*out_gS$SE),2)
	
#Upper Bound
out_gS$U=round(exp(log(out_gS$GMR)+(qt(1-out_gS$alpha1, DF))*out_gS$SE),2)


for (i in 1:nrow(df)) {
	if (out_gS[i,8] <= 0.8 | out_gS[i,9] >= 1.25) {	
	out_gS[i,10] = "STAGE 2"
	} else {
	out_gS[i,10] = "PASS EoStg1"
	}
}


### Stage 2
#Geometric Mean of Y (Not LN transformed Data) GMR K2
N2 = N1
N_Tot_gs = N1+N2
for (i in 1:nrow(df)) {
	if(out_gS[i,10] == "STAGE 2" ){
	
	out_gS[i,11] = N2		#N2	
	out_gS[i,12] = N1+N2	#Ntot
	
	out_r$GMR2_SD=exp(rowMeans(log(dfr[,1:N_Tot_gs]))) 
	out_t$GMR2_SD=exp(rowMeans(log(dft[,1:N_Tot_gs]))) 
	
	out_gS[i,13]=out_t[i,5]/out_r[i,5] 	#GMRk2
	
	out_r[i,6]=var(log(as.numeric(dfr[i,1:N_Tot_gs]))) 		#Vpoolk2 Test 
	out_t[i,6]=var(log(as.numeric(dft[i,1:N_Tot_gs]))) 		#Vpoolk2 Reference
	out_gS[i,14]=(out_t[i,6]+out_r[i,6])/2		#Vpoolk2
	out_gS[i,15]=sqrt(2*out_gS[i,14]/out_gS[i,12]/2) 				#SEk2

	#alpha2=0.0294 and DF
	out_gS[i,16] = 0.0294
	DF=out_gS[i,12]-2
		
	#t_apha_df2
	out_gS[i,17] = qt(1-out_gS[i,16], DF)

	out_gS[i,18]=round(exp(log(out_gS[i,13])-(qt(1-out_gS[i,16], DF))*out_gS[i,15]),2)	#Lower Bound
	out_gS[i,19]=round(exp(log(out_gS[i,13])+(qt(1-out_gS[i,16], DF))*out_gS[i,15]),2)	 #Upper Bound
	
	} else {
	out_gS[i,12] = N1	#Ntot
	}
}

#"EoStg2"
for (i in 1:nrow(df)) {
	if (out_gS[i,10] == "STAGE 2" & (out_gS[i,18] < 0.8 | out_gS[i,19] > 1.25)) {
	out_gS[i,20] = as.character("FAIL")
	} else if ( out_gS[i,10] == "STAGE 2" & out_gS[i,18] >= 0.8 & out_gS[i,19] <= 1.25) {
	out_gS[i,20] = as.character("PASS EoStg2")
	} else {
	out_gS[i,20] = out_gS[i,10]
	}
}

#OverAll
for (i in 1:nrow(df)) {
	if (out_gS[i,20] == "PASS EoStg2" |  out_gS[i,20] == "PASS EoStg1") {
	out_gS[i,21] = as.character("PASS")
	} else {
	out_gS[i,21] = as.character("FAIL")
	}
}

#Export
output_gS = paste0("C:\\Users\\georgi.dezhov\\Documents\\Dipl\\OUTPUT\\SD_", "_", N_arm, ".csv")
write.csv(out_gS,output_gS, row.names = TRUE)



tabl = c(min(out_gS[,12]), quantile(out_gS[,12],  probs = 0.05), median(out_gS[,12]),
		mean(out_gS[,12]),quantile(out_gS[,12],  probs = 0.95), max(out_gS[,12]), length(out_gS[,12][out_gS[,12]<=54]))
names(tabl)=c("Min", "5%", "Median", "Mean", "95%", "Max", "<=54")
tabl

SDtabl <- matrix(0, ncol = 7, 1)
SDtabl = data.frame(SDtabl)
colnames(SDtabl)=c("Min", "5%", "Median", "Mean", "95%", "Max", "<=54")
rownames(SDtabl) = "Total"
SDtabl[1,] = tabl
SDtabl
write.csv(SDtabl,"C:\\Users\\georgi.dezhov\\Documents\\Dipl\\OUTPUT\\SDtabl.csv", row.names = TRUE)
prop.table(out_gS[,10])
sd_plot = table(out_gS[,10])
names(sd_plot)[1] = "STAGE 1"

#####################################
Total_sample_Adaptive = out[,11]
# JPEG device
jpeg("Total_sample_SD.jpeg", quality = 75)

# Code
barplot(sd_plot)

# Close device
dev.off()
#####################################

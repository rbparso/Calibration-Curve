## Import Data
library(readxl)
data = read_excel('C:\\Users\\rbparso\\OneDrive - Emory University\\Emory\\Datasets\\Summary-Calibration-Curve-Emory&Tricore-FDA.xlsx', sheet="DataProcess")

data = as.data.frame(data)
names(data) = c("Dilution", "BD.N1", "BD.N2", "TRICORE.Cobas.ORF1", "TRICORE.Cobas.E",
   "Emory.CDC.N2", "Emory.CDC.SC2", "Emory.Cobas.ORF1", "Emory.Cobas.E", "ddPCR",
   "Quest.Cobas.ORF1a", "Quest.Cobas.E", "ClearDx.ORF1", "ClearDx.E", "Log.ddPCR2", "N2",
   "Cepheid.E", "Cepheid.N2", "Cepheid.RdRP", "Roche.6800.ORF1ab", "Roche.6800.E", "Log.ddPCR3", "N2.CT")
data$Log.ddPCR = log(data[,10])

## Runs for Roche E2
library(mcr)
x=log(unlist(data[,10]))
y=50-unlist(data[,5])
pb.reg = mcreg(x, y, method.reg = "PaBa")
pb.reg@para
MCResult.plot(pb.reg, x.lab="log(ddPCR)", y.lab=expression(paste("50 - Roche E2 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 7.87 Slope = 1.13",
    add.cor=FALSE, add.grid=FALSE)

data2 = data[-c(1:20,114),]
x2=log(unlist(data2[,10]))
y2=50-unlist(data2[,5])
pb.reg2 = mcreg(x2, y2, method.reg = "PaBa")
pb.reg2@para
dev.new()
MCResult.plot(pb.reg2, x.lab="log(ddPCR)", y.lab=expression(paste("50 - Roche E2 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outlier 1 and 2, Intercept = 8.55 Slope = 1.07",
    add.cor=FALSE, add.grid=FALSE)

cvx = cvy = rep(0,11)
mx = my = rep(0,11)
vx = vy = rep(0,11)
sx = sy = rep(0,11)
for (i in 1:11) {
   rows = which(data2[,1] == i+1)
   mx[i] = mean(data2[rows,24])
   vx[i] = var(data2[rows,24])
   sx[i] = sd(data2[rows,24])
   cvx[i] = sx[i]/mx[i]*100
   my[i] = mean(data2[rows,5])
   vy[i] = var(data2[rows,5])
   sy[i] = sd(data2[rows,5])
   cvy[i] = sy[i]/my[i]*100
}
delta = vy/vx
cvrat = cvy/cvx
sdrat = sy/sx

vx = c(round(vx,2), round(mean(vx),2))
vy = c(round(vy,2), round(mean(vy),2))
delta = c(round(delta,2), round(mean(delta),2))
sdrat = c(round(sdrat,9), round(mean(sdrat),2))

dilution = c(2:12, "Average")
roche_mat = matrix(c(dilution, vx, vy, delta, sdrat), ncol=5)
roche_dat = as.data.frame(roche_mat)
colnames(roche_dat) = c("Dilution", "Var X: Var(Log(ddPCR))", 
    "Var Y: Var(50 - (Roche E2 CT))", "Delta: Var(Y)/Var(X)", "CV(Y)/CV(X)")
library(r2rtf)
library(dplyr)
roche_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//roche_dat.rtf")

dem.reg = mcreg(x2, y2, method.reg = "Deming", error.ratio=(1/delta[12]))
dem.reg@para

sxx = cov(x2,x2)
syy = cov(y2,y2)
sxy = cov(x2,y2)
md = delta[12]
dem.slope1 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int1 = mean(y2) - dem.slope1*mean(x2)

dev.new()
MCResult.plot(dem.reg, x.lab="log(ddPCR)", y.lab=expression(paste("50- Roche E2 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 5.31, Intercept = 8.56 Slope = 1.06",
    add.cor=FALSE, add.grid=FALSE)

target = exp((20-dem.reg@para[1,1])/dem.reg@para[2,1])
target


## Runs for NAMAS data ##

xpb1=data[,22]
ypb1=50-data[,20]
pb.reg1 = mcreg(xpb1, ypb1, method.reg = "PaBa")
pb.reg1@para
dev.new()
MCResult.plot(pb.reg1, x.lab="log(ddPCR)", y.lab="50 - (Roche 6800 ORF1ab CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 7.56 Slope = 1.30",
    add.cor=FALSE, add.grid=FALSE)

data_orfab = data[-c(233,235,240,227,230,16,18),]
xpb11=data_orfab[,22]
ypb11=50-data_orfab[,20]
pb.reg11 = mcreg(xpb11, ypb11, method.reg = "PaBa")
pb.reg11@para
dev.new()
MCResult.plot(pb.reg11, x.lab="log(ddPCR)", y.lab="50 - (Roche 6800 ORF1ab CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outliers, Intercept = 7.49 Slope = 1.31",
    add.cor=FALSE, add.grid=FALSE)

vorfab=rep(0,12)
vx1=rep(0,12)
for (i in 1:12) {
   rows = which(data_orfab[,1] == i)
   vx1[i] = var(data_orfab[rows,22])
   vorfab[i] = var(data_orfab[rows,20], na.rm=TRUE)
}
deltaorfab = vorfab/vx1
deltaorfab = c(deltaorfab, mean(deltaorfab))

xd1 = data_orfab[,22]
yd1 = 50-data_orfab[,20]
dem.reg1 = mcreg(xd1, yd1, method.reg = "Deming", error.ratio=(1/deltaorfab[13]))
deltaorfab[13]
dem.reg1@para
dev.new()
MCResult.plot(dem.reg1, x.lab="log(ddPCR)", y.lab="50-(Roche 6800 ORF1ab CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 3.25, Intercept = 7.71 Slope = 1.29",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd1,xd1)
syy = cov(yd1,yd1)
sxy = cov(xd1,yd1)
md = deltaorfab[13]
dem.slope1 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int1 = mean(yd1) - dem.slope1*mean(xd1)

target1 = 50 - (dem.reg1@para[1,1] + dem.reg1@para[2,1]*log(target))
target1

xpb2=data[,22]
ypb2=50-data[,21]
pb.reg2 = mcreg(xpb2, ypb2, method.reg = "PaBa")
pb.reg2@para
dev.new()
MCResult.plot(pb.reg2, x.lab="log(ddPCR)", y.lab="50 - (Roche 6800 E CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 5.62 Slope = 1.39",
    add.cor=FALSE, add.grid=FALSE)

data_cobase = data[-c(233,235,240,227,230,16,18),]
xpb22=data_cobase[,22]
ypb22=50-data_cobase[,21]
pb.reg22 = mcreg(xpb22, ypb22, method.reg = "PaBa")
pb.reg22@para
dev.new()
MCResult.plot(pb.reg22, x.lab="log(ddPCR)", y.lab="50 - (Roche 6800 E CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outliers, Intercept = 5.56 Slope = 1.40",
    add.cor=FALSE, add.grid=FALSE)

vcobase=rep(0,12)
vx2=rep(0,12)
for (i in 1:12) {
   rows = which(data_cobase[,1] == i)
   vx2[i] = var(data_cobase[rows,22])
   vcobase[i] = var(data_cobase[rows,21], na.rm=TRUE)
}
deltacobase = vcobase/vx2
deltacobase = c(deltacobase, mean(deltacobase))

dilution = c(1:12, "Average")
namas_mat = matrix(c(dilution, deltaorfab, deltacobase), ncol=3)
namas_dat = as.data.frame(namas_mat)
colnames(roche_dat) = c("Dilution", "Delta Roche 6800 ORF1ab: Var(Y1)/Var(X)", 
    "Delta Roche 6800 E Gene: Var(Y2)/Var(X)")
library(r2rtf)
library(dplyr)
namas_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//namas_dat.rtf")

xd2 = data_cobase[,22]
yd2 = 50-data_cobase[,21]
dem.reg2 = mcreg(xd2, yd2, method.reg = "Deming", error.ratio=(1/deltacobase[13]))
deltacobase[13]
dem.reg2@para
dev.new()
MCResult.plot(dem.reg2, x.lab="log(ddPCR)", y.lab="50-(Roche 6800 Cobas E CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 3.89, Intercept = 5.74 Slope = 1.39",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd2,xd2)
syy = cov(yd2,yd2)
sxy = cov(xd2,yd2)
md = deltacobase[13]
dem.slope2 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int2 = mean(yd2) - dem.slope2*mean(xd2)

target2 = 50 - (dem.reg2@para[1,1] + dem.reg2@para[2,1]*log(target))
target2
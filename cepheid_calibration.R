## Import Data
library(readxl)
data = read_excel('C:\\Users\\rbparso\\OneDrive - Emory University\\Emory\\Datasets\\Summary-Calibration-Curve-Emory&Tricore-FDA.xlsx', sheet="DataProcess")

data = as.data.frame(data)
names(data) = c("Dilution", "BD.N1", "BD.N2", "TRICORE.Cobas.ORF1", "TRICORE.Cobas.E",
   "Emory.CDC.N2", "Emory.CDC.SC2", "Emory.Cobas.ORF1", "Emory.Cobas.E", "ddPCR",
   "Quest.Cobas.ORF1a", "Quest.Cobas.E", "ClearDx.ORF1", "ClearDx.E", "Log.ddPCR2", "N2",
   "Cepheid.E", "Cepheid.N2", "Cepheid.RdRP")
data$Log.ddPCR = log(data[,10])

## Calculate calibrated ddPCR using old reads ##
library(mcr)
x=log(unlist(data[,10]))
y=50-unlist(data[,5])
pb.reg = mcreg(x, y, method.reg = "PaBa")
pb.reg@para
MCResult.plot(pb.reg, x.lab="log(ddPCR)", y.lab="50 - Roche E2 Ct", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 7.87 Slope = 1.13",
    add.cor=FALSE, add.grid=FALSE)

data2 = data[-c(1:20,114),]
x2=log(unlist(data2[,10]))
y2=50-unlist(data2[,5])
pb.reg2 = mcreg(x2, y2, method.reg = "PaBa")
pb.reg2@para
dev.new()
MCResult.plot(pb.reg2, x.lab="log(ddPCR)", y.lab="50 - Roche E2 Ct", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outlier 1 and 2, Intercept = 8.55 Slope = 1.07",
    add.cor=FALSE, add.grid=FALSE)

cvx = cvy = rep(0,11)
mx = my = rep(0,11)
vx = vy = rep(0,11)
sx = sy = rep(0,11)
for (i in 1:11) {
   rows = which(data2[,1] == i+1)
   mx[i] = mean(data2[rows,20])
   vx[i] = var(data2[rows,20])
   sx[i] = sd(data2[rows,20])
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
MCResult.plot(dem.reg, x.lab="log(ddPCR)", y.lab="50- Roche E2 Ct", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 5.31, Intercept = 8.56 Slope = 1.06",
    add.cor=FALSE, add.grid=FALSE)

target = exp((20-dem.reg@para[1,1])/dem.reg@para[2,1])
target

## Runs for Cepheid data ##

xpb1=data[,20]
ypb1=50-data[,17]
pb.reg1 = mcreg(xpb1, ypb1, method.reg = "PaBa")
pb.reg1@para
dev.new()
MCResult.plot(pb.reg1, x.lab="log(ddPCR)", y.lab="50 - (Cepheid E CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 5.81 Slope = 1.52",
    add.cor=FALSE, add.grid=FALSE)

data_cephe = data[-c(230,227,222,208,146,110),]
xpb11=data_cephe[,20]
ypb11=50-data_cephe[,17]
pb.reg11 = mcreg(xpb11, ypb11, method.reg = "PaBa")
pb.reg11@para
dev.new()
MCResult.plot(pb.reg11, x.lab="log(ddPCR)", y.lab="50 - (Cepheid E CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outliers, Intercept = 5.77 Slope = 1.52",
    add.cor=FALSE, add.grid=FALSE)

vcephe=rep(0,12)
vx1=rep(0,12)
for (i in 1:12) {
   rows = which(data_cephe[,1] == i)
   vx1[i] = var(data_cephe[rows,20])
   vcephe[i] = var(data_cephe[rows,17], na.rm=TRUE)
}
deltacephe = vcephe/vx1
deltacephe = c(deltacephe, mean(deltacephe))

xd1 = data_cephe[,20]
yd1 = 50-data_cephe[,17]
dem.reg1 = mcreg(xd1, yd1, method.reg = "Deming", error.ratio=(1/deltacephe[13]))
deltacephe[13]
dem.reg1@para
dev.new()
MCResult.plot(dem.reg1, x.lab="log(ddPCR)", y.lab="50-(Cepheid E CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 7.51, Intercept = 5.92 Slope = 1.51",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd1,xd1)
syy = cov(yd1,yd1)
sxy = cov(xd1,yd1)
md = deltacephe[13]
dem.slope1 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int1 = mean(yd1) - dem.slope1*mean(xd1)

target1 = 50 - (dem.reg1@para[1,1] + dem.reg1@para[2,1]*log(target))
target1

xpb2=data[,20]
ypb2=50-data[,18]
pb.reg2 = mcreg(xpb2, ypb2, method.reg = "PaBa")
pb.reg2@para
dev.new()
MCResult.plot(pb.reg2, x.lab="log(ddPCR)", y.lab="50 - (Cepheid N2 CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 2.13 Slope = 1.51",
    add.cor=FALSE, add.grid=FALSE)

data_cephn = data[-c(230,227,224,222,146,110),]
xpb22=data_cephn[,20]
ypb22=50-data_cephn[,18]
pb.reg22 = mcreg(xpb22, ypb22, method.reg = "PaBa")
pb.reg22@para
dev.new()
MCResult.plot(pb.reg22, x.lab="log(ddPCR)", y.lab="50 - (Cepheid N2 CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outliers, Intercept = 2.07 Slope = 1.52",
    add.cor=FALSE, add.grid=FALSE)

vcephn=rep(0,12)
vx2=rep(0,12)
for (i in 1:12) {
   rows = which(data_cephn[,1] == i)
   vx2[i] = var(data_cephn[rows,20])
   vcephn[i] = var(data_cephn[rows,18], na.rm=TRUE)
}
deltacephn = vcephn/vx2
deltacephn = c(deltacephn, mean(deltacephn))

xd2 = data_cephn[,20]
yd2 = 50-data_cephn[,18]
dem.reg2 = mcreg(xd2, yd2, method.reg = "Deming", error.ratio=(1/deltacephn[13]))
deltacephn[13]
dem.reg2@para
dev.new()
MCResult.plot(dem.reg2, x.lab="log(ddPCR)", y.lab="50-(Cepheid N2 CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 9.49, Intercept = 2.18 Slope = 1.51",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd2,xd2)
syy = cov(yd2,yd2)
sxy = cov(xd2,yd2)
md = deltacephn[13]
dem.slope2 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int2 = mean(yd2) - dem.slope2*mean(xd2)

target2 = 50 - (dem.reg2@para[1,1] + dem.reg2@para[2,1]*log(target))
target2

xpb3=data[,20]
ypb3=50-data[,19]
pb.reg3 = mcreg(xpb3, ypb3, method.reg = "PaBa")
pb.reg3@para
dev.new()
MCResult.plot(pb.reg3, x.lab="log(ddPCR)", y.lab="50 - (Cepheid RdRP CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 3.04 Slope = 1.54",
    add.cor=FALSE, add.grid=FALSE)

data_cephr = data[-c(230,208,146,110),]
xpb33=data_cephr[,20]
ypb33=50-data_cephr[,19]
pb.reg33 = mcreg(xpb33, ypb33, method.reg = "PaBa")
pb.reg33@para
dev.new()
MCResult.plot(pb.reg33, x.lab="log(ddPCR)", y.lab="50 - (Cepheid RdRP CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outliers, Intercept = 3.06 Slope = 1.54",
    add.cor=FALSE, add.grid=FALSE)

vcephr=rep(0,12)
vx3=rep(0,12)
for (i in 1:12) {
   rows = which(data_cephr[,1] == i)
   vx3[i] = var(data_cephr[rows,20])
   vcephr[i] = var(data_cephr[rows,19], na.rm=TRUE)
}
deltacephr = vcephr/vx3
deltacephr = c(deltacephr, mean(deltacephr))

dilution = c(1:12, "Average")
ceph_mat = matrix(c(dilution, deltacephe, deltacephn, deltacephr), ncol=4)
ceph_dat = as.data.frame(ceph_mat)
colnames(ceph_dat) = c("Dilution",
    "Cepheid E Delta: Var(Y1)/Var(X)", "Cepheid N2 Delta: Var(Y2)/Var(X)", 
    "Cepheid RdRP Delta: Var(Y2)/Var(X)")
ceph_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//ceph_dat.rtf")


xd3 = data_cephr[,20]
yd3 = 50-data_cephr[,19]
dem.reg3 = mcreg(xd3, yd3, method.reg = "Deming", error.ratio=(1/deltacephr[13]))
dem.reg3@para
deltacephr[13]
dev.new()
MCResult.plot(dem.reg3, x.lab="log(ddPCR)", y.lab="50-(Cepheid RdRP CT)", points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 8.20, Intercept = 3.28 Slope = 1.52",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd3,xd3)
syy = cov(yd3,yd3)
sxy = cov(xd3,yd3)
md = deltacephr[13]
dem.slope3 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int3 = mean(yd3) - dem.slope3*mean(xd3)

target3 = 50 - (dem.reg3@para[1,1] + dem.reg3@para[2,1]*log(target))
target3
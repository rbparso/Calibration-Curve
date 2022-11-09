## Import Data
library(readxl)
data = read_excel('C:\\Users\\rbparso\\OneDrive - Emory University\\Emory\\Datasets\\Summary-Calibration-Curve-Emory&Tricore-FDA.xlsx', sheet="DataProcess")

data = as.data.frame(data)
names(data) = c("Dilution", "BD.N1", "BD.N2", "TRICORE.Cobas.ORF1", "TRICORE.Cobas.E",
   "Emory.CDC.N2", "Emory.CDC.SC2", "Emory.Cobas.ORF1", "Emory.Cobas.E", "ddPCR",
   "Quest.Cobas.ORF1a", "Quest.Cobas.E", "ClearDx.ORF1", "ClearDx.E", "Log.ddPCR2", "N2",
   "Cepheid E", "Cepheid N2", "Cepheid RdRP")
data$Log.ddPCR = log(data[,10])

## Runs for Roche E2
library(mcr)
x=log(unlist(data[,10]))
y=50-unlist(data[,5])
pb.reg = mcreg(x, y, method.reg = "PaBa")
pb.reg@para
MCResult.plot(pb.reg, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - Roche cobas\U00AE E ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 7.87 Slope = 1.13",
    add.cor=FALSE, add.grid=FALSE)

data2 = data[-c(1:20,114),]
x2=log(unlist(data2[,10]))
y2=50-unlist(data2[,5])
pb.reg2 = mcreg(x2, y2, method.reg = "PaBa")
pb.reg2@para
dev.new()
MCResult.plot(pb.reg2, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - Roche cobas\U00AE E ",C[T],sep="")), points.pch=16,
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
sdrat = c(round(sdrat,2), round(mean(sdrat),2))

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
MCResult.plot(dem.reg, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50- Roche cobas\U00AE E ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 5.31, Intercept = 8.56 Slope = 1.06",
    add.cor=FALSE, add.grid=FALSE)

target = exp((20-dem.reg@para[1,1])/dem.reg@para[2,1])
target

## Runs for BD N1/N2
data3 = data[-c(139,229),]
y3=50-data3[,2]
x3=data3[,24]
pb.reg3 = mcreg(x3, y3, method.reg = "PaBa")
pb.reg3@para
dev.new()
MCResult.plot(pb.reg3, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - BD N1 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 8.42 Slope = 1.36",
    add.cor=FALSE, add.grid=FALSE)

data4 = data[-c(139,223,229),]
y4=50-data4[,3]
x4=data4[,24]
pb.reg4 = mcreg(x4, y4, method.reg = "PaBa")
pb.reg4@para
dev.new()
MCResult.plot(pb.reg4, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - BD N2 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 7.46 Slope = 1.38",
    add.cor=FALSE, add.grid=FALSE)

data5 = data[-c(139,229,223,224,230),]
x5 = data5[,24]
y5 = 50-data5[,2]
pb.reg5 = mcreg(x5, y5, method.reg = "PaBa")
pb.reg5@para
dev.new()
MCResult.plot(pb.reg5, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - BD N1 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outlier, Intercept = 8.35 Slope = 1.37",
    add.cor=FALSE, add.grid=FALSE)

data6 = data[-c(139,229,223,230),]
x6 = data6[,24]
y6 = 50-data6[,3]
pb.reg6 = mcreg(x6, y6, method.reg = "PaBa")
pb.reg6@para
dev.new()
MCResult.plot(pb.reg6, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - BD N2 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outlier, Intercept = 7.43 Slope = 1.38",
    add.cor=FALSE, add.grid=FALSE)

mn1=rep(0,12)
vn1=rep(0,12)
sn1=rep(0,12)
mx2=rep(0,12)
vx2=rep(0,12)
sx2=rep(0,12)
for (i in 1:12) {
   rows = which(data5[,1] == i)
   mx2[i] = mean(data5[rows,24], na.rm=TRUE)
   vx2[i] = var(data5[rows,24], na.rm=TRUE)
   sx2[i] = sd(data5[rows,24], na.rm=TRUE)
   mn1[i] = mean(data5[rows,2], na.rm=TRUE)
   vn1[i] = var(data5[rows,2], na.rm=TRUE)
   sn1[i] = sd(data5[rows,2], na.rm=TRUE)
}
deltan1 = vn1/vx2
sdratn1 = sn1/sx2

mn2=rep(0,12)
vn2=rep(0,12)
sn2=rep(0,12)
mx2=rep(0,12)
vx2=rep(0,12)
sx2=rep(0,12)
for (i in 1:12) {
   rows = which(data5[,1] == i)
   mx2[i] = mean(data5[rows,24])
   vx2[i] = var(data5[rows,24])
   sx2[i] = sd(data5[rows,24])
   mn2[i] = mean(data5[rows,3], na.rm=TRUE)
   vn2[i] = var(data5[rows,3], na.rm=TRUE)
   sn2[i] = sd(data5[rows,3], na.rm=TRUE)
}
deltan2 = vn2/vx2
sdratn2 = sn2/sx2

deltan1 = c(round(deltan1,3), round(mean(deltan1),3))
deltan2 = c(round(deltan2,3), round(mean(deltan2),3))
sdratn1 = c(round(sdratn1,3), round(mean(sdratn1),3))
sdratn2 = c(round(sdratn2,3), round(mean(sdratn2),3))

dilution = c(1:12, "Average")
bd_mat = matrix(c(dilution, deltan1, deltan2, sdratn1, sdratn2), ncol=5)
bd_dat = as.data.frame(bd_mat)
colnames(bd_dat) = c("Dilution", "BD Max N1 Delta: Var(Y1)/Var(X)", 
    "BD Max N2 Delta: Var(Y2)/Var(X)", "BD Max N1 CV(Y1)/CV(X)", 
    "BD Max N2 CV(Y2)/CV(X)")
library(r2rtf)
library(dplyr)
bd_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//bd_dat.rtf")

xd2 = data5[,24]
yd2 = 50-data5[,2]
dem.reg2 = mcreg(xd2, yd2, method.reg = "Deming", error.ratio=(1/mean(deltan1[1:12])))
dem.reg2@para
dev.new()
MCResult.plot(dem.reg2, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50-(BD N1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 9.406, Intercept = 8.41 Slope = 1.37",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd2,xd2)
syy = cov(yd2,yd2)
sxy = cov(xd2,yd2)
md = deltan1[13]
dem.slope2 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int2 = mean(yd2) - dem.slope2*mean(xd2)

xd3 = data5[,24]
yd3 = 50-data5[,3]
dem.reg3 = mcreg(xd3, yd3, method.reg = "Deming", error.ratio=(1/mean(deltan2[1:12])))
dem.reg3@para
dev.new()
MCResult.plot(dem.reg3, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50-(BD N2 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 19.246, Intercept = 7.55 Slope = 1.38",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd3,xd3)
syy = cov(yd3,yd3)
sxy = cov(xd3,yd3)
md = deltan2[13]
dem.slope3 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int3 = mean(yd3) - dem.slope3*mean(xd3)

target2 = 50 - (dem.reg2@para[1,1] + dem.reg2@para[2,1]*log(target))
target2
target3 = 50 - (dem.reg3@para[1,1] + dem.reg3@para[2,1]*log(target))
target3

## Runs for Emory CDC
x7=data[,24]
y7=50-data[,6]
pb.reg7 = mcreg(x7, y7, method.reg = "PaBa")
pb.reg7@para
dev.new()
MCResult.plot(pb.reg7, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - CDC N2 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 5.57 Slope = 1.45",
    add.cor=FALSE, add.grid=FALSE)

x8=data[,24]
y8=50-data[,7]
pb.reg8 = mcreg(x8, y8, method.reg = "PaBa")
pb.reg8@para
dev.new()
MCResult.plot(pb.reg8, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - CDC SC2 ",C[T],sep="")),points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 8.08 Slope = 1.21",
    add.cor=FALSE, add.grid=FALSE)

data9 = data[-c(222,224,227,229,230),]
x9=data9[,24]
y9=50-data9[,7]
pb.reg9 = mcreg(x9, y9, method.reg = "PaBa")
pb.reg9@para
dev.new()
MCResult.plot(pb.reg9, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50 - CDC SC2 ",C[T],sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 7.98 Slope = 1.22",
    add.cor=FALSE, add.grid=FALSE)

mcdc=rep(0,12)
vcdc=rep(0,12)
scdc=rep(0,12)
mx3=rep(0,12)
vx3=rep(0,12)
sx3=rep(0,12)
for (i in 1:12) {
   rows = which(data9[,1] == i)
   mx3[i] = mean(data9[rows,24])
   vx3[i] = var(data9[rows,24])
   sx3[i] = sd(data9[rows,24])
   mcdc[i] = mean(data9[rows,6], na.rm=TRUE)
   vcdc[i] = var(data9[rows,6], na.rm=TRUE)
   scdc[i] = sd(data9[rows,6], na.rm=TRUE)
}

msc2=rep(0,12)
vsc2=rep(0,12)
ssc2=rep(0,12)
mx4=rep(0,12)
vx4=rep(0,12)
sx4=rep(0,12)
for (i in 1:12) {
   rows = which(data9[,1] == i)
   mx4[i] = mean(data9[rows,24])
   vx4[i] = var(data9[rows,24])
   sx4[i] = sd(data9[rows,24])
   msc2[i] = mean(data9[rows,7], na.rm=TRUE)
   vsc2[i] = var(data9[rows,7], na.rm=TRUE)
   ssc2[i] = sd(data9[rows,7], na.rm=TRUE)
}
deltacdc = vcdc/vx3
sdratcdc = scdc/sx3
deltasc2 = vsc2/vx4
sdratsc2 = ssc2/sx4

vx4 = c(round(vx4,2), round(mean(vx4),2))
deltacdc = c(round(deltacdc,2), round(mean(deltacdc),2))
deltasc2 = c(round(deltasc2,2), round(mean(deltasc2),2))
sdratcdc = c(round(sdratcdc,2), round(mean(sdratcdc),2))
sdratsc2 = c(round(sdratsc2,2), round(mean(sdratsc2),2))

dilution = c(1:12, "Average")
em_mat = matrix(c(dilution, vx4, deltacdc, deltasc2, sdratcdc, sdratsc2), ncol=6)
em_dat = as.data.frame(em_mat)
colnames(em_dat) = c("Dilution", "Var(X): Var(Log(ddPCR))",
    "CDC N2 Delta: Var(Y1)/Var(X)", "CDC SC2 Delta: Var(Y2)/Var(X)", 
    "CV(Y1)/CV(X)", "CV(Y2)/CV(X)")
em_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//em_dat.rtf")

xd4 = data[,24]
yd4 = 50-data[,6]
dem.reg4 = mcreg(xd4, yd4, method.reg = "Deming", error.ratio=(1/mean(deltacdc[1:12])))
dem.reg4@para
dev.new()
MCResult.plot(dem.reg4, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50-(CDC N2 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 3.11, Intercept = 5.68 Slope = 1.44",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd4,xd4)
syy = cov(yd4,yd4)
sxy = cov(xd4,yd4)
md = deltacdc[13]
dem.slope4 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int4 = mean(yd4) - dem.slope4*mean(xd4)

xd5 = data9[,24]
yd5 = 50-data9[,7]
dem.reg5 = mcreg(xd5, yd5, method.reg = "Deming", error.ratio=(1/mean(deltasc2[1:12])))
dem.reg5@para
dev.new()
MCResult.plot(dem.reg5, x.lab="Natural-Log Transformed ddPCR", 
    y.lab=expression(paste("50-(CDC SC2 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 1.97, Intercept = 8.28 Slope = 1.2",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd5,xd5)
syy = cov(yd5,yd5)
sxy = cov(xd5,yd5)
md = deltasc2[13]
dem.slope5 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int5 = mean(yd5) - dem.slope5*mean(xd5)

target4 = 50 - (dem.reg4@para[1,1] + dem.reg4@para[2,1]*log(target))
target4
target5 = 50 - (dem.reg5@para[1,1] + dem.reg5@para[2,1]*log(target))
target5

## Runs for TriCore/Emory ORF1 ##

xpb1=data[,24]
ypb1=50-data[,4]
pb.reg8 = mcreg(xpb1, ypb1, method.reg = "PaBa")
pb.reg8@para
dev.new()
MCResult.plot(pb.reg8, x.lab="log(ddPCR)", y.lab=expression(paste("50 - (TriCore ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 6.14 Slope = 1.31",
    add.cor=FALSE, add.grid=FALSE)

vtorf=rep(0,12)
vx5=rep(0,12)
for (i in 1:12) {
   rows = which(data[,1] == i)
   vx5[i] = var(data[rows,24])
   vtorf[i] = var(data[rows,4], na.rm=TRUE)
}
deltatorf = vtorf/vx5
deltatorf = c(deltatorf, mean(deltatorf))

xd6 = data[,24]
yd6 = 50-data[,4]
dem.reg6 = mcreg(xd6, yd6, method.reg = "Deming", error.ratio=(1/deltatorf[13]))
deltatorf[13]
dem.reg6@para
dev.new()
MCResult.plot(dem.reg6, x.lab="log(ddPCR)", y.lab=expression(paste("50-(TriCore ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 6.30, Intercept = 5.97 Slope = 1.32",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd6,xd6)
syy = cov(yd6,yd6)
sxy = cov(xd6,yd6)
md = deltatorf[13]
dem.slope6 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int6 = mean(yd6) - dem.slope6*mean(xd6)

target6 = 50 - (dem.reg6@para[1,1] + dem.reg6@para[2,1]*log(target))
target6

xpb2=data[,20]
ypb2=50-data[,8]
pb.reg9 = mcreg(xpb2, ypb2, method.reg = "PaBa")
pb.reg9@para
dev.new()
MCResult.plot(pb.reg9, x.lab="log(ddPCR)", y.lab=expression(paste("50 - (Emory ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 6.82 Slope = 1.38",
    add.cor=FALSE, add.grid=FALSE)

veorf=rep(0,12)
vx6=rep(0,12)
for (i in 1:12) {
   rows = which(data[,1] == i)
   vx6[i] = var(data[rows,20])
   veorf[i] = var(data[rows,8], na.rm=TRUE)
}
deltaeorf = veorf/vx6
deltaeorf = c(deltaeorf, mean(deltaeorf))

xd7 = data[,20]
yd7 = 50-data[,8]
dem.reg7 = mcreg(xd7, yd7, method.reg = "Deming", error.ratio=(1/deltaeorf[13]))
deltaeorf[13]
dem.reg7@para
dev.new()
MCResult.plot(dem.reg7, x.lab="log(ddPCR)", y.lab=expression(paste("50-(Emory ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 4.82, Intercept = 7.17 Slope = 1.35",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd7,xd7)
syy = cov(yd7,yd7)
sxy = cov(xd7,yd7)
md = deltaeorf[13]
dem.slope7 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int7 = mean(yd7) - dem.slope7*mean(xd7)

target7 = 50 - (dem.reg7@para[1,1] + dem.reg7@para[2,1]*log(target))
target7

## Run for Emory Cobas E ##

xpb3=data[,20]
ypb3=50-data[,9]
pb.reg10 = mcreg(xpb3, ypb3, method.reg = "PaBa")
pb.reg10@para
dev.new()
MCResult.plot(pb.reg10, x.lab="log(ddPCR)", y.lab=expression(paste("50 - (Emory Cobas E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 5.66 Slope = 1.40",
    add.cor=FALSE, add.grid=FALSE)

vece=rep(0,12)
vx7=rep(0,12)
for (i in 1:12) {
   rows = which(data[,1] == i)
   vx7[i] = var(data[rows,20])
   vece[i] = var(data[rows,9], na.rm=TRUE)
}
deltaece = vece/vx7
deltaece = c(deltaece, mean(deltaece))

dilution = c(1:12, "Average")
triem_mat = matrix(c(dilution, deltatorf, deltaeorf, deltaece), ncol=4)
triem_dat = as.data.frame(triem_mat)
colnames(triem_dat) = c("Dilution",
    "TriCore ORF1 Delta: Var(Y1)/Var(X)", "Emory ORF1 Delta: Var(Y2)/Var(X)", 
    "Emory Cobas E Delta: Var(Y2)/Var(X)")
triem_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//triem_dat.rtf")


xd8 = data[,20]
yd8 = 50-data[,9]
dem.reg8 = mcreg(xd8, yd8, method.reg = "Deming", error.ratio=(1/deltaece[13]))
dem.reg8@para
deltaece[13]
dev.new()
MCResult.plot(dem.reg8, x.lab="log(ddPCR)", y.lab=expression(paste("50-(Emory Cobas E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 8.89, Intercept = 5.96 Slope = 1.37",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd8,xd8)
syy = cov(yd8,yd8)
sxy = cov(xd8,yd8)
md = deltaece[13]
dem.slope8 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int8 = mean(yd8) - dem.slope8*mean(xd8)

target8 = 50 - (dem.reg8@para[1,1] + dem.reg8@para[2,1]*log(target))
target8

## Run for Quest ORF1 ##

data10 = data[complete.cases(data[,11]),]
xpb4=data10[,20]
ypb4=50-data10[,11]
pb.reg11 = mcreg(xpb4, ypb4, method.reg = "PaBa")
pb.reg11@para
dev.new()
MCResult.plot(pb.reg11, x.lab="log(ddPCR)", y.lab=expression(paste("50 - (Quest ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 6.54 Slope = 1.30",
    add.cor=FALSE, add.grid=FALSE)

vqorf=rep(0,12)
vx8=rep(0,12)
for (i in 1:12) {
   rows = which(data[,1] == i)
   vx8[i] = var(data[rows,20])
   vqorf[i] = var(data[rows,11], na.rm=TRUE)
}
deltaqorf = vqorf/vx8
deltaqorf = c(deltaqorf, mean(deltaqorf))

xd9 = data10[,20]
yd9 = 50-data10[,11]
dem.reg9 = mcreg(xd9, yd9, method.reg = "Deming", error.ratio=(1/deltaqorf[13]))
deltaqorf[13]
dem.reg9@para
dev.new()
MCResult.plot(dem.reg9, x.lab="log(ddPCR)", y.lab=expression(paste("50-(Quest ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 27.33, Intercept = 7.00 Slope = 1.26",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd9,xd9)
syy = cov(yd9,yd9)
sxy = cov(xd9,yd9)
md = deltaqorf[13]
dem.slope9 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int9 = mean(yd9) - dem.slope9*mean(xd9)

target9 = 50 - (dem.reg9@para[1,1] + dem.reg9@para[2,1]*log(target))
target9

## Run for Quest Cobas E ##

dataqe = data[complete.cases(data[,12]),]
xpb5=dataqe[,20]
ypb5=50-dataqe[,12]
pb.reg12 = mcreg(xpb5, ypb5, method.reg = "PaBa")
pb.reg12@para
dev.new()
MCResult.plot(pb.reg12, x.lab="log(ddPCR)", y.lab=expression(paste("50 - (Quest Cobas E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 4.27 Slope = 1.41",
    add.cor=FALSE, add.grid=FALSE)

vqce=rep(0,12)
vx9=rep(0,12)
for (i in 1:12) {
   rows = which(dataqe[,1] == i)
   vx9[i] = var(dataqe[rows,20])
   vqce[i] = var(dataqe[rows,12], na.rm=TRUE)
}
deltaqce = vqce/vx9
deltaqce = c(deltaqce, mean(deltaqce))

dilution = c(1:12, "Average")
quest_mat = matrix(c(dilution, deltaqorf, deltaqce), ncol=3)
quest_dat = as.data.frame(quest_mat)
colnames(quest_dat) = c("Dilution",
    "Quest ORF1a Delta: Var(Y1)/Var(X)", "Quest Cobas E Delta: Var(Y2)/Var(X)")
quest_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//quest_dat.rtf")

xd10 = dataqe[,20]
yd10 = 50-dataqe[,12]
dem.reg10 = mcreg(xd10, yd10, method.reg = "Deming", error.ratio=(1/deltaqce[13]))
deltaqce[13]
dem.reg10@para
dev.new()
MCResult.plot(dem.reg10, x.lab="log(ddPCR)", y.lab=expression(paste("50-(Quest Cobas E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 30.40, Intercept = 4.71 Slope = 1.37",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd10,xd10)
syy = cov(yd10,yd10)
sxy = cov(xd10,yd10)
md = deltaqce[13]
dem.slope10 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int10 = mean(yd10) - dem.slope10*mean(xd10)

target10 = 50 - (dem.reg10@para[1,1] + dem.reg10@para[2,1]*log(target))
target10

## Run for ClearDx ORF1a ##

xpb6=data[,20]
ypb6=50-data[,13]
pb.reg13 = mcreg(xpb6, ypb6, method.reg = "PaBa")
pb.reg13@para
dev.new()
MCResult.plot(pb.reg13, x.lab="log(ddPCR)", y.lab=expression(paste("50 - (ClearDx ORF1a ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 5.85 Slope = 1.40",
    add.cor=FALSE, add.grid=FALSE)

vcorf=rep(0,12)
vx10=rep(0,12)
for (i in 1:12) {
   rows = which(data[,1] == i)
   vx10[i] = var(data[rows,20])
   vcorf[i] = var(data[rows,13], na.rm=TRUE)
}
deltacorf = vcorf/vx10
deltacorf = c(deltacorf, mean(deltacorf))

xd11 = data[,20]
yd11 = 50-data[,13]
dem.reg11 = mcreg(xd11, yd11, method.reg = "Deming", error.ratio=(1/deltacorf[13]))
deltacorf[13]
dem.reg11@para
dev.new()
MCResult.plot(dem.reg11, x.lab="log(ddPCR)", y.lab=expression(paste("50-(ClearDx ORF1a ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 3.14, Intercept = 5.99 Slope = 1.38",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd11,xd11)
syy = cov(yd11,yd11)
sxy = cov(xd11,yd11)
md = deltacorf[13]
dem.slope11 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int11 = mean(yd11) - dem.slope11*mean(xd11)

target11 = 50 - (dem.reg11@para[1,1] + dem.reg11@para[2,1]*log(target))
target11

## Run for ClearDx E ##

xpb7=data[,20]
ypb7=50-data[,14]
pb.reg14 = mcreg(xpb7, ypb7, method.reg = "PaBa")
pb.reg14@para
dev.new()
MCResult.plot(pb.reg14, x.lab="log(ddPCR)", y.lab=expression(paste("50 - (ClearDx E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 7.31 Slope = 1.26",
    add.cor=FALSE, add.grid=FALSE)

vce=rep(0,12)
vx11=rep(0,12)
for (i in 1:12) {
   rows = which(data[,1] == i)
   vx11[i] = var(data[rows,20])
   vce[i] = var(data[rows,14], na.rm=TRUE)
}
deltace = vce/vx11
deltace = c(deltace, mean(deltace))

dilution = c(1:12, "Average")
clear_mat = matrix(c(dilution, deltacorf, deltace), ncol=3)
clear_dat = as.data.frame(clear_mat)
colnames(clear_dat) = c("Dilution",
    "ClearDx ORF1a Delta: Var(Y1)/Var(X)", "ClearDx E Delta: Var(Y2)/Var(X)")
clear_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//clear_dat.rtf")

xd12 = data[,20]
yd12 = 50-data[,14]
dem.reg12 = mcreg(xd12, yd12, method.reg = "Deming", error.ratio=(1/deltace[13]))
deltace[13]
dem.reg12@para
dev.new()
MCResult.plot(dem.reg12, x.lab="log(ddPCR)", y.lab=expression(paste("50-(ClearDx E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 2.61, Intercept = 7.34 Slope = 1.26",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd12,xd12)
syy = cov(yd12,yd12)
sxy = cov(xd12,yd12)
md = deltace[13]
dem.slope12 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int12 = mean(yd12) - dem.slope12*mean(xd12)

target12 = 50 - (dem.reg12@para[1,1] + dem.reg12@para[2,1]*log(target))
target12

## ClearDx ORF1 using new ddPCR ##

xpb10=data[,15]
ypb10=50-data[,13]
pb.reg17 = mcreg(xpb10, ypb10, method.reg = "PaBa")
pb.reg17@para
dev.new()
MCResult.plot(pb.reg17, x.lab="log(ddPCR2)", y.lab=expression(paste("50 - (ClearDx ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 6.38 Slope = 1.38",
    add.cor=FALSE, add.grid=FALSE)

data_coag = data[-c(235,233,240,230,227,16,18),]
xpb11=data_coag[,15]
ypb11=50-data_coag[,13]
pb.reg18 = mcreg(xpb11, ypb11, method.reg = "PaBa")
pb.reg18@para
dev.new()
MCResult.plot(pb.reg18, x.lab="log(ddPCR2)", y.lab=expression(paste("50 - (ClearDx ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outliers, Intercept = 6.30 Slope = 1.38",
    add.cor=FALSE, add.grid=FALSE)

vcoag=rep(0,12)
vx15=rep(0,12)
for (i in 1:12) {
   rows = which(data_coag[,1] == i)
   vx15[i] = var(data_coag[rows,15])
   vcoag[i] = var(data_coag[rows,12], na.rm=TRUE)
}
deltacoag = vcoag/vx15
deltacoag = c(deltacoag, mean(deltacoag))

xd15 = data_coag[,15]
yd15 = 50-data_coag[,13]
dem.reg15 = mcreg(xd15, yd15, method.reg = "Deming", error.ratio=(1/deltacoag[13]))
deltacoag[13]
dem.reg15@para
dev.new()
MCResult.plot(dem.reg15, x.lab="log(ddPCR2)", y.lab=expression(paste("50-(ClearDx ORF1 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 33.75, Intercept = 6.28 Slope = 1.39",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd15,xd15)
syy = cov(yd15,yd15)
sxy = cov(xd15,yd15)
md = deltacoag[13]
dem.slope15 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int15 = mean(yd15) - dem.slope15*mean(xd15)

target17 = 50 - (dem.reg15@para[1,1] + dem.reg15@para[2,1]*log(target))
target17

## ClearDx E Ag using new ddPCR ##

xpb8=data[,15]
ypb8=50-data[,14]
pb.reg15 = mcreg(xpb8, ypb8, method.reg = "PaBa")
pb.reg15@para
dev.new()
MCResult.plot(pb.reg15, x.lab="log(ddPCR2)", y.lab=expression(paste("50 - (ClearDx E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 7.84 Slope = 1.23",
    add.cor=FALSE, add.grid=FALSE)

data_ceag = data[-c(235,233,240,230,227,16,18),]
xpb8=data_ceag[,15]
ypb8=50-data_ceag[,14]
pb.reg15 = mcreg(xpb8, ypb8, method.reg = "PaBa")
pb.reg15@para
dev.new()
MCResult.plot(pb.reg15, x.lab="log(ddPCR2)", y.lab=expression(paste("50 - (ClearDx E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outliers, Intercept = 7.73 Slope = 1.24",
    add.cor=FALSE, add.grid=FALSE)

vceag=rep(0,12)
vx12=rep(0,12)
for (i in 1:12) {
   rows = which(data_ceag[,1] == i)
   vx12[i] = var(data_ceag[rows,15])
   vceag[i] = var(data_ceag[rows,14], na.rm=TRUE)
}
deltaceag = vceag/vx12
deltaceag = c(deltaceag, mean(deltaceag))

xd13 = data_ceag[,15]
yd13 = 50-data_ceag[,14]
dem.reg13 = mcreg(xd13, yd13, method.reg = "Deming", error.ratio=(1/deltaceag[13]))
deltaceag[13]
dem.reg13@para
dev.new()
MCResult.plot(dem.reg13, x.lab="log(ddPCR2)", y.lab=expression(paste("50-(ClearDx E ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 3.17, Intercept = 7.56 Slope = 1.27",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd13,xd13)
syy = cov(yd13,yd13)
sxy = cov(xd13,yd13)
md = deltaceag[13]
dem.slope13 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int13 = mean(yd13) - dem.slope13*mean(xd13)

target15 = 50 - (dem.reg13@para[1,1] + dem.reg13@para[2,1]*log(target))
target15

target13 = exp((20-dem.reg13@para[1,1])/dem.reg13@para[2,1])
target13

## Emory N2 Ag using new ddPCR ##

xpb9=data[,15]
ypb9=50-data[,16]
pb.reg16 = mcreg(xpb9, ypb9, method.reg = "PaBa")
pb.reg16@para
dev.new()
MCResult.plot(pb.reg16, x.lab="log(ddPCR2)", y.lab=expression(paste("50 - (Emory N2 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Intercept = 5.75 Slope = 1.45",
    add.cor=FALSE, add.grid=FALSE)

data_n2ag = data[c(-235,-233,-230,-240,-227,-16, -18),]
xpb9=data_n2ag[,15]
ypb9=50-data_n2ag[,16]
pb.reg16 = mcreg(xpb9, ypb9, method.reg = "PaBa")
pb.reg16@para
dev.new()
MCResult.plot(pb.reg16, x.lab="log(ddPCR2)", y.lab=expression(paste("50 - (Emory N2 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.border=TRUE, ci.border.col=c("red", "black"),
    main="Passing BabLok Regression, Removed Outliers, Intercept = 5.71 Slope = 1.45",
    add.cor=FALSE, add.grid=FALSE)

vn2ag=rep(0,12)
vx13=rep(0,12)
for (i in 1:12) {
   rows = which(data_n2ag[,1] == i)
   vx13[i] = var(data_n2ag[rows,15])
   vn2ag[i] = var(data_n2ag[rows,16], na.rm=TRUE)
}
deltan2ag = vn2ag/vx13
deltan2ag = c(deltan2ag, mean(deltan2ag))

dilution = c(1:12, "Average")
newag_mat = matrix(c(dilution, deltacoag, deltaceag, deltan2ag), ncol=4)
newag_dat = as.data.frame(newag_mat)
colnames(newag_dat) = c("Dilution", "ClearDx ORF1 Delta: Var(Y1)/Var(X)",
    "ClearDx E Delta: Var(Y2)/Var(X)", "Emory N2 Delta: Var(Y3)/Var(X)")
newag_dat %>%
  rtf_body() %>%
  rtf_encode() %>%
  write_rtf(file="C://Users//rbparso//OneDrive - Emory University//Emory//Output//ct//newag_dat.rtf")


xd14 = data_n2ag[,15]
yd14 = 50-data_n2ag[,16]
dem.reg14 = mcreg(xd14, yd14, method.reg = "Deming", error.ratio=(1/deltan2ag[13]))
deltan2ag[13]
dem.reg14@para
dev.new()
MCResult.plot(dem.reg14, x.lab="log(ddPCR2)", y.lab=expression(paste("50-(N2 ",C[T],")",sep="")), points.pch=16,
    reg.col="darkgray", identity=FALSE, ci.area=FALSE, 
    main="Deming Reg, Delta = 2.59, Intercept = 5.77 Slope = 1.44",
    add.cor=FALSE, add.grid=FALSE)

sxx = cov(xd14,xd14)
syy = cov(yd14,yd14)
sxy = cov(xd14,yd14)
md = deltan2ag[13]
dem.slope14 = (syy-md*sxx+sqrt((syy-md*sxx)^2+4*md*sxy^2))/(2*sxy)
dem.int14 = mean(yd14) - dem.slope14*mean(xd14)

target14 = exp((20-dem.reg14@para[1,1])/dem.reg14@para[2,1])
target14

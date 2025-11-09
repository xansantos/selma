library(data.table)
library(foreach)
library(ggplot2)
library(doParallel)
library(MuMIn)


#length dependent catch comparison curves for each species (Standardbootstrapping)
#wichtig: diese 3 Dateien müssen in R eingeladen sein, damit Skript laufen kann

source("/home/santos/Documents/4_PROJECTS/STELLA2/Hannah/2024/glmAVG/F-glmAvg.R")
source("/home/santos/Documents/4_PROJECTS/STELLA2/Hannah/2024/glmAVG/b-glmAvg.R")
source("/home/santos/Documents/4_PROJECTS/STELLA2/Hannah/2024/glmAVG/B-glmAvg-doparallel.R")



# Data processing ---------------------------------------------------------


dat_raw <- fread("C:\\Users\\santos\\Nextcloud\\software_transfer\\selma_beta_1\\catch_stat.csv") 


summary(dat_raw)

length_columns<-colnames(dat_raw)[grep( "L_",colnames(dat_raw))]

#Haul, Perle, Art bleiben als Zeilen, Länge muss zur Zeile gemacht werden
dat_long<-melt(data=dat_raw,id.vars=c("Haul","String","Perle","Art"), measure.vars=length_columns,
     variable.name = "Length", value.name = "n",na.rm = FALSE, variable.factor = FALSE)

#Reference: Standardnetz; Test: Perlennetz
#wenn Perle =0 dann Referenz; sonst Test
dat_long[, Net:=ifelse(Perle==0, "Reference", "Test")]

dat_long[, Length:=type.convert(substr(Length, 3, 4), as.is=T)]

dat_wide<-dcast(dat_long, Haul+String+Art+Length~Net, fun.aggregate = sum,value.var ="n")

#Umbenennung der Spaltennamen in Haul, sp, l, n2 & n1
setnames(dat_wide, c("haul","string","sp","l","n_2","n_1"))

#q1 und q2=1, da die Länge von allen Fischen gemessen wurde
dat_wide[,':='(n=n_1+n_2, q_1=1, q_2=1)]

#überall, wo Länge = 0 beim Standard- und Perlennetz ist, wird gelöscht
dat_wide<-dat_wide[n>0]


# Catch comparison PLE ----------------------------------------------------

X_ple<-dat_wide[sp=="PLE",.(haul,string,l,n_1,n_2,q_1,q_2)]

X_ple<-setDF(X_ple[order(haul,string,l)])

pooled_ple<-pooldata(X=X_ple)


# to run the glms, first source these scripts -----------------------------


## bootstrap model, confidence intervals calculated directly on response values
mod_ple<- selma_boot(DT_in = X_ple,
                                         lengthRange=c(19,45),
                                         typePredict="response",
                                         ranking="AICc",
                                         sampler="bootsamplerGN",
                                         exportfuns=c("pooldata","selma"),
                                         ncores=4,
                                         B=1000)

#original ggplot2 plots

p_cc_ple<-ggplot(data=mod_ple[[3]],aes(x=l))+ ylim(0,1)+
  geom_hline(yintercept = 0.5,col="darkgreen")+
  geom_line( aes(y=mu),col=2,lwd=1.5) + xlab("fish length (cm)") +ylab("catch propotion in pearlnet")+
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.2)+
  geom_point(data=mod_ple[[2]],aes(x=l,y=p_l,size=n_l), alpha=.5,col="darkblue")+
  #ggtitle("average catch comparison curves of plaice between a standard-gillnet and a pearlnet")+
  ggtitle("PLE Analysis")+
  theme(plot.title = element_text(color="red", size=14, face="bold.italic"))

print(p_cc_ple)

#length dependent catch ratio curve

p_cr_ple<-ggplot(data=mod_ple[[4]],aes(x=l))+
  geom_hline(yintercept = 1.0,col="darkgreen")+
  #xlim(15,50)+
  coord_cartesian(xlim=c(19,46),ylim=c(0,3))+
  geom_line( aes(y=mu),col=2,lwd=1.5) + xlab("fish length (cm)") +ylab("Catch ratio of TEST to CONTROL ")+ theme(text = element_text(size=rel(4)))+
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.2)

print(p_cr_ple)

# bespoken figures  ----------------------------------------------------
library(data.table)


### some colors:
col_p_up <- colors()[24]  

polygonN<-rgb(.5, .5, .5,alpha=.9)

polygonCC<-rgb(.2, .2, .2,alpha=.2)

col.points1<-rgb(.8, .8, .8,alpha=.4)

### to automatically specify labels sequence for length frequency y axis
f_yax<-function(yst){ 
  stp=if(yst<10) 2 else {
    if(yst>=10 & yst<26) 5 else{
      if(yst>=26 & yst<75) 10  else{
        if(yst>=75 & yst<100) 25 else{
          if(yst>=100 & yst<150) 25 else{
            if(yst>=150 & yst<300) 50  else{
              if(yst>=300 & yst<501) 100  else{ 250
              }}}}}}}
  
  return(stp)}
par(mar=c(1, 1, 5, 6),
    mfrow=c(1,1),
    oma = c(4, 4, 0, 1),
    cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

plot(mu~l,data=mod_ple[[3]],type="n",
     main="",lwd=2,yaxt="n",xaxt="n",lty=1,bty="n",xlim=c(15,45),ylim=c(0,1200),ylab="",xlab="")
axis(1, at = seq(20, 75, by = 5),line=0,las=1)
axis(4, at = seq(0, 1200, by = f_yax(1200)),line=0,las=2)
polygon(c(rev(pooled_ple$l), pooled_ple$l), 
        c(rev(pooled_ple$n_2), rep(0,length(pooled_ple$n_2))),
        col=polygonN, 
        border = NA)

lines(n_1~l,data=pooled_ple,type="l",lwd=2,lty=2)

par(new=TRUE)

plot(mu~l,data=mod_ple[[3]],type = 'n',bty="n",ylim=c(0,1),xlim=c(15,45),ylab="",xlab="Length",main="",las=1)
polygon(c(rev(mod_ple[[3]][["l"]]), mod_ple[[3]][["l"]]),
        c(rev(mod_ple[[3]][["uci"]]), mod_ple[[3]][["lci"]]), 
        col = polygonCC, border = NA)
axis(2,at=seq(0,1,0.1),las=2)
lines(mod_ple[[3]][["mu"]]~ mod_ple[[3]][["l"]],lty=1,lwd=2,col=1)
lines(mod_ple[[3]][["lci"]]~ mod_ple[[3]][["l"]],lty = 1,lwd=2,col="white")
lines(mod_ple[[3]][["uci"]]~ mod_ple[[3]][["l"]],lty = 1,lwd=2,col="white")
abline(h=0.5,col="darkgreen",lty=3,lwd=2)
abline(v=25,col="darkgreen",lty=3,lwd=2)
points(p_l~l,data=mod_ple[[2]], pch=21, col="red", bg=1)




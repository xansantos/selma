library(data.table)
library(foreach)
library(ggplot2)
library(doParallel)
library(MuMIn)

lapply(grep("funs",list.files()), function(i) source(list.files()[i],echo=TRUE))

COD<-fread("./data/COD.csv")



# to run the glms, first source these scripts -----------------------------


## bootstrap model, confidence intervals calculated directly on response values
mod_cod<- selma_boot(DT_in = COD,
                             lengthRange=c(19,45),
                             typePredict="response",
                             ranking="AICc",
                             sampler="bootsamplerGN",
                             exportfuns=c("pooldata","selma"),
                             ncores=4,
                             B=100
                     )

mod_ple<- selma_boot(DT_in = PLE,
                     lengthRange=c(19,45),
                     typePredict="response",
                     ranking="AICc",
                     sampler="bootsamplerGN",
                     exportfuns=c("pooldata","selma"),
                     ncores=4,
                     B=1000
)

saveRDS(mod_ple,file="./data/mod_ple.rds")

#original ggplot2 plots

p_cc_cod<-ggplot(data=mod_cod[[3]],aes(x=l))+ ylim(0,1)+
  geom_hline(yintercept = 0.5,col="darkgreen")+
  geom_line( aes(y=mu),col=2,lwd=1.5) + xlab("fish length (cm)") +ylab("catch propotion in pearlnet")+
  geom_ribbon(aes(ymin=lci, ymax=uci), alpha=0.2)+
  geom_point(data=mod_cod[[2]],aes(x=l,y=p_l,size=n_l), alpha=.5,col="darkblue")+
  #ggtitle("average catch comparison curves of plaice between a standard-gillnet and a pearlnet")+
  ggtitle("COD Analysis")+
  theme(plot.title = element_text(color="red", size=14, face="bold.italic"))

print(p_cc_cod)

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




#0. Opcje i biblioteki ------------------------------------------------------

options(numeric=50,scipen = 50);library(ggplot2);library(data.table);library(plotly);
library(gridExtra);library(zoo);library(scales);library(grid)
##
#1. Przykładowe trajkektorie analityczne ------------------------------------
N<-2^15
T=1
ileSymulacji=3

zbior<-do.call(rbind,lapply(1:ileSymulacji,function(x){
br<-c(0,rnorm(N-1,0,sqrt(T/N)))
  data.frame(
    Xanalityczne=exp(-(-cumsum(br)+sqrt(2))^2),
    t=seq(0,T,length=N),
    symulacja=x
  )
}))

ggplot(zbior)+
  geom_line(aes(t,Xanalityczne,colour=factor(symulacja)))+
  theme_bw()
#2. Przykładowa trajektoria eulera -----------------------------------------

N<-2^15
T<-3
br<-c(0,rnorm(N,0,sqrt(T/N)))
Xeuler<-Xmillstein<-c()
Xeuler[1]<-Xmillstein[1]<-exp(-2)


for(x in 2:N){
  if(!(Xeuler[length(Xeuler)]<1&Xeuler[length(Xeuler)]<1&
       Xmillstein[length(Xmillstein)]<1&Xmillstein[length(Xmillstein)]>0)){break}
   if(!(Xeuler[length(Xeuler)]<=0|Xeuler[length(Xeuler)]>=1)){
  Xeuler[x]<-Xeuler[x-1]+
    -Xeuler[x-1]*(2*log(Xeuler[x-1])+1)*T/N+
    2*(Xeuler[x-1])*sqrt(-log(Xeuler[x-1]))*br[x]}
  if(!(Xmillstein[length(Xmillstein)]<=0|Xmillstein[length(Xmillstein)]>=1)){
  Xmillstein[x]<-Xmillstein[x-1]+
    -Xmillstein[x-1]*(2*log(Xmillstein[x-1])+1)*
    T/N+2*(Xmillstein[x-1])*sqrt(-log(Xmillstein[x-1]))*br[x]+
    Xmillstein[x-1]*(-2*log(Xmillstein[x-1])-1)*
    (br[x]^2-T/N)}
}

br5<-c(0,rollapply(br[-1],by=2^5,width=2^5,FUN=sum))
Xeuler.5<-Xmillstein.5<-c()
Xeuler.5[1]<-Xmillstein.5[1]<-exp(-2)

for(x in 2:(N/2^5+1)){
  if(!(Xeuler.5[length(Xeuler.5)]<1&Xeuler.5[length(Xeuler.5)]<1&
       Xmillstein.5[length(Xmillstein.5)]<1&Xmillstein.5[length(Xmillstein.5)]>0)){break}
  if(!(Xeuler.5[length(Xeuler.5)]<=0|Xeuler.5[length(Xeuler.5)]>=1)){
    Xeuler.5[x]<-Xeuler.5[x-1]+
      -Xeuler.5[x-1]*(2*log(Xeuler.5[x-1])+1)*4*2^5*T/N+
      2*(Xeuler.5[x-1])*sqrt(-log(Xeuler.5[x-1]))*br5[x]}
  if(!(Xmillstein.5[length(Xmillstein.5)]<=0|Xmillstein.5[length(Xmillstein.5)]>=1)){
    Xmillstein.5[x]<-Xmillstein.5[x-1]+
      -Xmillstein.5[x-1]*(2*log(Xmillstein.5[x-1])+1)*
      4*2^5*T/N+2*(Xmillstein.5[x-1])*sqrt(-log(Xmillstein.5[x-1]))*br5[x]+
      Xmillstein.5[x-1]*(-2*log(Xmillstein.5[x-1])-1)*
      (br5[x]^2-4*2^5*T/N)}
}

br10<-c(0,rollapply(br[-1],by=2^10,width=2^10,FUN=sum))
Xeuler.10<-Xmillstein.10<-c()
Xeuler.10[1]<-Xmillstein.10[1]<-exp(-2)

for(x in 2:(N/2^10+1)){
  if(!(Xeuler.10[length(Xeuler.10)]<1&Xeuler.10[length(Xeuler.10)]<1&
       Xmillstein.10[length(Xmillstein.10)]<1&Xmillstein.10[length(Xmillstein.10)]>0)){break}
  if(!(Xeuler.10[length(Xeuler.10)]<=0|Xeuler.10[length(Xeuler.10)]>=1)){
    Xeuler.10[x]<-Xeuler.10[x-1]+
      -Xeuler.10[x-1]*(2*log(Xeuler.10[x-1])+1)*2^10*T/N+
      2*(Xeuler.10[x-1])*sqrt(-log(Xeuler.10[x-1]))*br10[x]}
  if(!(Xmillstein.10[length(Xmillstein.10)]<=0|Xmillstein.10[length(Xmillstein.10)]>=1)){
    Xmillstein.10[x]<-Xmillstein.10[x-1]+
      -Xmillstein.10[x-1]*(2*log(Xmillstein.10[x-1])+1)*
      2^10*T/N+2*(Xmillstein.10[x-1])*sqrt(-log(Xmillstein.10[x-1]))*br10[x]+
      Xmillstein.10[x-1]*(-2*log(Xmillstein.10[x-1])-1)*
      (br10[x]^2-2^10*T/N)}
}


Xanalityczne<-exp(-(-cumsum(br)+sqrt(2))^2)
zbior<-data.frame(
  Proces=c(Xeuler,Xanalityczne,Xmillstein),
  Rodzaj=c(rep('X_t euler',length(Xeuler)),
           rep('X_t analityczne',length(Xanalityczne)),
           rep('X_t millstein',length(Xmillstein))),
  t=c(seq(0,T,length=N+1)[1:length(Xeuler)],
      seq(0,T,length=N+1)[1:length(Xanalityczne)],
      seq(0,T,length=N+1)[1:length(Xmillstein)])
)

zbior2<-data.frame(
  Proces=c(Xanalityczne,Xeuler,Xeuler.5,Xeuler.10),
  Rodzaj=c(rep('X_t analityczne',length(Xanalityczne)),
           rep('X_t euler, 2^15',length(Xeuler)),
           rep('X_t euler, 2^10',length(Xeuler.5)),
           rep('X_t euler, 2^5',length(Xeuler.10))),
  t=c(seq(0,T,length=N+1)[1:length(Xanalityczne)],
      seq(0,T,length=N+1)[1:length(Xeuler)],
      seq(0,T,length=N/2^5+1)[1:length(Xeuler.5)],
      seq(0,T,length=N/2^10+1)[1:length(Xeuler.10)])
)

ggplot()+
  geom_line(data=zbior,aes(t,Proces,colour=factor(Rodzaj)),alpha=.8)+
  scale_colour_manual(values = rainbow(3))+
  geom_line(data=data.frame(x=seq(0,T,length=N)[length(Xeuler)],y=c(0,1)),aes(x=x,y=y),
            linetype=2,colour=rainbow(3)[2])+
  geom_line(data=data.frame(x=seq(0,T,length=N)[length(Xmillstein)],y=c(0,1)),aes(x=x,y=y),
            linetype=2,colour=rainbow(3)[3])+
  theme_bw()+
  theme(text = element_text(size=20))+
  labs(colour='Rodzaj')

kolory<-hue_pal()(3)

ggplot()+
  geom_line(data=zbior2,aes(t,Proces,colour=factor(Rodzaj)),lwd=1)+
  scale_colour_manual(values = c('black',kolory))+
  geom_line(data=data.frame(x=seq(0,T,length=N+1)[length(Xeuler)],y=c(0,1)),aes(x=x,y=y),
            linetype=2,colour=kolory[2],lwd=1)+
  geom_line(data=data.frame(x=seq(0,T,length=N/2^5+1)[length(Xeuler.5)],y=c(0,1)),aes(x=x,y=y),
            linetype=2,colour=kolory[1],lwd=1)+
  geom_line(data=data.frame(x=seq(0,T,length=N/2^10+1)[length(Xeuler.10)],y=c(0,1)),aes(x=x,y=y),
            linetype=2,colour=kolory[3],lwd=1)+
  theme_bw()+
  theme(text = element_text(size=20))+
  labs(colour='Rodzaj')
  ##
#3. Millstein ---------------------------------------------------------------
N<-2^15
T<-3
br<-c(0,rnorm(N-1,0,sqrt(T/N)))
X<-c()
X[1]<-exp(-2)


for(x in 2:N){
  if(X[x-1]<=0|X[x-1]>=1){break}
  X[x]<-X[x-1]+
    -X[x-1]*(2*log(X[x-1])+1)*T/N+2*(X[x-1])*sqrt(-log(X[x-1]))*br[x]+
    (X[x-1])*sqrt(-log(X[x-1]))*(-2*sqrt(-log(X[x-1]))+1/sqrt(-log(X[x-1])))*(br[x]^2-T/N)
}

Xanalityczne<-exp(-(-cumsum(br)+sqrt(2))^2)
zbior<-data.frame(
  Proces=c(X,Xanalityczne),
  Rodzaj=c(rep('X',length(X)),rep('X_t analityczne',length(br))),
  t=c(seq(0,T,length=N)[1:length(X)],seq(0,T,length=N)[1:length(br)])
)

ggplot(zbior)+
  geom_line(aes(t,Proces,colour=factor(Rodzaj,levels =rev(unique(Rodzaj)))),alpha=.8)+
  scale_colour_manual(values=c('#330000','red'))+
  labs(colour='rodzaj')+
  theme_bw()
##
#4. funkcja momentów --------------------------------------------------------

EX<-function(t){
  sqrt(1/(1+2*t))*exp(-2/(1+2*t))
}
M2<-function(t){
  sqrt(1/(1+4*t))*exp(-4/(1+4*t))
}
VarX<-function(t){M2(t)-EX(t)^2}
EX3<- function(t) {
  sqrt(1/(1+6*t))*exp(-6/(1+6*t))
}
##
#5. deszcz trajektori -------------------------------------------------------


N<-2^8
T<-200
zbior<-do.call(rbind,lapply(1:400,function(x){
    br<-c(0,rnorm(N,0,sqrt(T/N)))
  return(data.frame(
    t=seq(0,T,length=N+1),
    symulacja=x,
    X=exp(-(-cumsum(br)+sqrt(2))^2)
  ))
}))
ggplot(zbior)+
  geom_line(aes(t,X,group=factor(symulacja)),alpha=.05,colour='#000033',lwd=2)+
  theme_bw()
##
#6. zestawienie z funkcją mt ------------------------------------------------

N<-2^12
T<-1
part<-100
Symulacji<-1000/part
zbior<-do.call(rbind,lapply(1:Symulacji,function(x){

             df<-do.call(rbind,lapply(1:part,function(y){
               br<-c(0,rnorm(N,0,sqrt(T/N)))
               Xeuler<-Xmillstein<-c()
               Xeuler[1]<-Xmillstein[1]<-exp(-2)
               
               
               for(x in 2:N){
                 if(!(Xeuler[length(Xeuler)]<1&Xeuler[length(Xeuler)]>0&
                      Xmillstein[length(Xmillstein)]<1&Xmillstein[length(Xmillstein)]>0)){break}
                  if(!(Xeuler[length(Xeuler)]<=0|Xeuler[length(Xeuler)]>=1)){
                   Xeuler[x]<-Xeuler[x-1]+
                     -Xeuler[x-1]*(2*log((Xeuler[x-1]))+1)*T/N+
                     2*(Xeuler[x-1])*sqrt(-log(Xeuler[x-1]))*br[x]}
                 if(!(Xmillstein[length(Xmillstein)]<=0|Xmillstein[length(Xmillstein)]>=1)){
                   Xmillstein[x]<-Xmillstein[x-1]+
                     -Xmillstein[x-1]*(2*log(Xmillstein[x-1])+1)*
                     T/N+2*(Xmillstein[x-1])*sqrt(-log(Xmillstein[x-1]))*br[x]+
                     Xmillstein[x-1]*(-2*log(Xmillstein[x-1])-1)*
                     (br[x]^2-T/N)}
               }
               Xanalityczne<-exp(-(-cumsum(br)[-(N+1)]+sqrt(2))^2)
               zbior<-data.table(
                 Proces=c(Xeuler,Xanalityczne,Xmillstein),
                 Rodzaj=c(rep('X_t euler',length(Xeuler)),
                          rep('X_t analityczne',length(Xanalityczne)),
                          rep('X_t millstein',length(Xmillstein))),
                 t=c(seq(0,T,length=N)[1:length(Xeuler)],
                     seq(0,T,length=N)[1:length(Xanalityczne)],
                     seq(0,T,length=N)[1:length(Xmillstein)])
               )
               return(zbior)
             }));  return(df[,list(
               EX=mean(Proces),
               var=mean(Proces^2)-mean(Proces)^2,
               Ex3=mean(Proces^3),
               waga=length(Proces)
             ),by=list(Rodzaj,t)])}))
            
zbior<-data.table(zbior)

zbiorMomenty<-zbior[,list(
  EX=weighted.mean(EX,waga),
  var=weighted.mean(var,waga),
  Ex3=weighted.mean(Ex3,waga)
),by=list(Rodzaj,t)]


a<-ggplot(zbiorMomenty)+
  geom_line(aes(t,EX,colour=Rodzaj),lwd=1)+
  stat_function(fun=EX,colour='black',lwd=1)+
  theme_bw()+theme(text=element_text(size=16))+
  guides(colour=F)

b<-ggplot(zbiorMomenty)+
  geom_line(aes(t,var,colour=Rodzaj),lwd=1)+
  stat_function(fun=VarX,colour='black',lwd=1)+
  theme_bw()+theme(text=element_text(size=16))+
  guides(colour=F)

c<-ggplot(zbiorMomenty)+
  geom_line(aes(t,Ex3,colour=Rodzaj),lwd=1)+
  stat_function(fun=EX3,colour='black',lwd=1)+
  theme_bw()+theme(text=element_text(size=16))+
  theme(legend.justification=c(.95,0),
        legend.position=c(.95,0),
        legend.background = element_rect(fill="gray90", 
                                         size=.5, 
                                         linetype="dotted"))
main=textGrob('N=2^12,symulacji:1000',gp=gpar(fontsize=20,font=3))
grid.arrange(a,b,c,nrow=1,top=main)

##
#7. Kalkukacja momentów i gęstości ------------------------------------------
  
N<-2^15

zbior2<-do.call(rbind,lapply(1:10000,function(x){
br<-c(0,rnorm(N,0,sqrt(10/N)))
Xanalityczne<-exp(-(-cumsum(br)+sqrt(2))^2)
return(matrix(Xanalityczne[floor(seq(1/10*N,N,by=1/10*N))],nrow=1))
}))
for(i in 1:10){
  assign(paste('gestosc',i),density(zbior2[,i])$y)
}

zbior4<-data.frame(
  X=sapply(1:10,function(x)rep(x,512))[1:5120],
  Y=rep(seq(0,1,length=512),10),
  Z=sapply(1:10,function(x){get(paste('gestosc',x))})[1:5120]
)
plot_ly(data=zbior4, x = ~X, y = ~Y, z = ~Z,split = ~X, type = "scatter3d", mode = "lines")

ggplot(zbior4)+
  geom_line(aes(Y,Z,colour=factor(X)))
##
#8. Sprawdzenie gęstości analitycznej ---------------------------------------

#zdefiniowanie funkcji gęstości
gest<-function(t,x){
  1/(2*x*sqrt(-log(x)*2*pi*t))*exp(-(sqrt(-log(x))+sqrt(2))^2*1/(2*t))*(exp(2*sqrt(-2*log(x))/t)+1)
}
#zdefiniowanie punktów w których liczymy gęstość(siatka na płaszczyźnie)
x<-rep(seq(.01,.99,by=.01),10)
t<-do.call('c',lapply(seq(.2,2,length=10)^2,function(x)rep(x,99)))
#zapisanie współrzędnych do rysowania w zbiorze
Gestosci<-data.frame(
  x=x,
  t=t,
  z=gest(t,x)
)
#rysowanie ewolucji gęstości w 3d
plot_ly(data=Gestosci, x = ~t, y = ~x, z = ~z,
        split = ~paste(sqrt(t),'^2'),
        type = "scatter3d",
        mode = "lines",line=list(width=5))
#rysowanie ewolucji gęstości zrzutowanych na X_t x f(X_t)
ggplot(Gestosci)+
  geom_line(aes(x,z,group=t,colour=factor(paste(sqrt(t),'^2'))))+
  labs(colour='t')+
  theme_bw()
##
#9. Badanie rozkładu --------------------------------------------------------

N<-2^12
symulacji=1000
T<-c(1,2,8)



for(t in T){
  assign(paste('t.rowne.',t,sep=''),
         do.call(rbind,lapply(1:symulacji,function(n){
           br<-c(0,rnorm(N,0,sqrt(t/N)))
           Xeuler<-Xmillstein<-c()
           Xeuler[1]<-Xmillstein[1]<-exp(-2)
    for(x in 2:N){
      if(!(Xeuler[length(Xeuler)]<1&Xeuler[length(Xeuler)]<1&
           Xmillstein[length(Xmillstein)]<1&Xmillstein[length(Xmillstein)]>0)){
        break}
      if(!(Xeuler[length(Xeuler)]<=0|Xeuler[length(Xeuler)]>=1)){
        Xeuler[x]<-Xeuler[x-1]+
          -Xeuler[x-1]*(2*log(Xeuler[x-1])+1)*t/N+
          2*(Xeuler[x-1])*sqrt(-log(Xeuler[x-1]))*br[x]}
      if(!(Xmillstein[length(Xmillstein)]<=0|Xmillstein[length(Xmillstein)]>=1)){
        Xmillstein[x]<-Xmillstein[x-1]+
          -Xmillstein[x-1]*(2*log(Xmillstein[x-1])+1)*
          t/N+2*(Xmillstein[x-1])*sqrt(-log(Xmillstein[x-1]))*br[x]+
          Xmillstein[x-1]*(-2*log(Xmillstein[x-1])-1)*
          (br[x]^2-t/N)}
    }
               return(data.frame(
                 wartosc=c(Xeuler[N],Xmillstein[N]),
                 aproksymacja=c('Euler','Millstein')
           ))})))}

df<-subset(t.rowne.8,aproksymacja=='Euler')
a<-ggplot(df,aes(wartosc))+
  geom_histogram(aes(y=..density..),
                 fill='#000033',
                 alpha=.8,
                 bins=12)+
  stat_function(fun=function(x)
    gest(8,x),xlim = c(.01,.99),colour='red',lwd=1)+
  theme_bw()+
  ggtitle(paste('t=8,',sum(!is.na(df$wartosc)),'trajektori'))+
  theme(text=element_text(size=18))+
  labs(x='X_t')

df<-subset(t.rowne.2,aproksymacja=='Euler')
b<-ggplot(df,aes(wartosc))+
  geom_histogram(aes(y=..density..),
                 fill='#000033',
                 alpha=.8,
                 bins=12)+
  stat_function(fun=function(x)
    gest(2,x),xlim = c(.01,.99),colour='red',lwd=1)+
  theme_bw()+
  ggtitle(paste('t=2,',sum(!is.na(df$wartosc)),'trajektori'))+
  theme(text=element_text(size=18))+
  labs(x='X_t')

df<-subset(t.rowne.1,aproksymacja=='Euler')
c<-ggplot(df,aes(wartosc))+
  geom_histogram(aes(y=..density..),
                 fill='#000033',
                 alpha=.8,
                 bins=12)+
  stat_function(fun=function(x)
    gest(1,x),xlim = c(.01,.99),colour='red',lwd=1)+
  theme_bw()+
  ggtitle(paste('t=1,',sum(!is.na(df$wartosc)),'trajektori'))+
  theme(text=element_text(size=18))+
  labs(x='X_t')

grid.arrange(c,b,a,nrow=1)
#10. Błąd aproksymacji -------------------------------------------------------
Ny<-2^(5:10)
symulacji=2000
T<-1

zbiorBledow<-do.call(rbind,lapply(Ny,function(N){
  do.call(rbind,lapply(1:symulacji,function(n){
    br<-c(0,rnorm(N,0,sqrt(T/N)))
    Xeuler<-Xmillstein<-c()
    Xeuler[1]<-Xmillstein[1]<-exp(-2)
    for(x in 2:N){
      if(!(Xeuler[length(Xeuler)]<1&Xeuler[length(Xeuler)]<1&
           Xmillstein[length(Xmillstein)]<1&Xmillstein[length(Xmillstein)]>0)){
        break}
      if(!(Xeuler[length(Xeuler)]<=0|Xeuler[length(Xeuler)]>=1)){
        Xeuler[x]<-Xeuler[x-1]+
          -Xeuler[x-1]*(2*log(Xeuler[x-1])+1)*T/N+
          2*(Xeuler[x-1])*sqrt(-log(Xeuler[x-1]))*br[x]}
      if(!(Xmillstein[length(Xmillstein)]<=0|Xmillstein[length(Xmillstein)]>=1)){
        Xmillstein[x]<-Xmillstein[x-1]+
          -Xmillstein[x-1]*(2*log(Xmillstein[x-1])+1)*
          T/N+2*(Xmillstein[x-1])*sqrt(-log(Xmillstein[x-1]))*br[x]+
          Xmillstein[x-1]*(-2*log(Xmillstein[x-1])-1)*
          (br[x]^2-T/N)}
    }
    X_t<-exp(-(sum(br)-sqrt(2))^2)
    return(data.frame(
      blad=c(abs(Xeuler[N]-X_t),abs(Xmillstein[N]-X_t)),
      aproksymacja=c('Euler','Millstein'),
      n=T/N
    ))}))}))

zbiorBledow<-aggregate(blad~aproksymacja+n,
                       data=zbiorBledow,
                       FUN=mean)

#zestawienie z f=1/2x
ggplot()+
  geom_point(data=zbiorBledow,
             aes(x=n,y=blad,colour=factor(aproksymacja)),size=3)+
  scale_x_log10("dt",
                breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x)))+
  scale_y_log10("log(błąd)",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  geom_line(data=data.frame(x=c(0.001,.1),y=c(0.0305,.08)),
            aes(x,y),size=1)+
  theme_bw()+
  theme(text=element_text(size=18))+
  labs(colour='Metoda')+
  ggtitle('Średni błąd aproksymacji, 5000 symulacji.\nSkala log-log')

#dodanie lini regresji

ggplot(data=zbiorBledow,
       aes(x=n,y=blad,colour=factor(aproksymacja)))+
  geom_point(size=3)+
  scale_x_log10("dt",
                breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x)))+
  scale_y_log10("log(błąd)",
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  geom_smooth(method = 'lm',se=F)+
  theme_bw()+
  theme(text=element_text(size=18))+
  labs(colour='Metoda')+
  ggtitle('Średni błąd aproksymacji, 5000 symulacji.\nSkala log-log')
#
lm(log(zbiorBledow$blad[zbiorBledow$aproksymacja=='Euler'])~
     log(zbiorBledow$n[zbiorBledow$aproksymacja=='Euler']))$coefficients[2]->aEuler

lm(log(zbiorBledow$blad[zbiorBledow$aproksymacja=='Millstein'])~
     log(zbiorBledow$n[zbiorBledow$aproksymacja=='Millstein']))$coefficients[2]->aMillstein
cat("gamma wyestymowane euler:",aEuler,'\n',"gamma wyestymowane millstein:",aMillstein,'\n')

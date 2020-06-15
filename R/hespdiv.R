spatial.analysis<-function(data,method=NA,variation=NA,metric=NA,criteria=NA,
                           C.cond=0,E.cond=0,N.cond=0,S.cond=0,divisions=NULL,lim=NULL,
                           knot.density.X=knot.density.X,knot.density.Y=knot.density.Y,curve.iterations,
                           correction.term=0.05,null.models=T,seed.t=round(runif(1,0,9999),0),test.n) {
  ##### pirmas data stulpelis turi buti rusis - faktorius, antras X, trecias Y., dependencies = dplyr
  #S.cond pateikti kaip proporcija ploto
  require(dplyr)
  require(pracma)
  library(gstat)
  library(sp)
  library(spatstat)
  #duomenys pirminiai
  {
    plot(data$X,data$Y,col=as.numeric(data$rusis))
    ch<-chull(data$X, data$Y)
    ids<-c(ch,ch[1])
    x<-data$X[ids]
    y<-data$Y[ids]
    lines(x,y)
    origins.chull<-data.frame(X=x,Y=y)
    lines(origins.chull)
    rims<-list(origins.chull)
    blokai<-data.frame(performance=0,iteration=0,saknis=0)

    split.reliability<-numeric()
    split.performance<-numeric()
    split.reliability2<-numeric()
    n.splits<-numeric()
    ave.split.abE<-numeric()
  }
  S.cond<-abs(polyarea(x,y))*S.cond
  splits<-numeric()
  checks<-list()
  original.quality<-numeric()
  iteration<-1
  #motinine rekursyvine funkcija
  spatial_div<-function(samp.dat,method,variation,
                        metric,criteria,C.cond=0,E.cond=0,N.cond=0,S.cond=0,
                        root=2,divisions=NULL,rims,knot.density.X=knot.density.X,
                        knot.density.Y=knot.density.Y,curve.iterations=curve.iterations,correction.term,split.reliability,
                        splits,split.performance,split.reliability2,checks,null.models=null.models,n.splits=n.splits,seed.t=seed.t,
                        ave.split.abE=ave.split.abE,test.n=test.n){
    #testuojamas plotas
    testid<-length(rims)
    margins<-rims[[testid]]
    original.qual<- -2*entropija(samp.dat[,1])
    iteration<<- iteration +1
    #grafikas pirminis nupaisomas
    {
      plot(data$X,data$Y,col=data$rusis)
      points(samp.dat$X,samp.dat$Y,pch=19)
      centras<-poly_center(margins[,1],margins[,1])
      points(centras[1],centras[2],col=3,pch=19)
      lines(rims[[1]])
    }
    print(c("testid: ", testid))
    # padalinimai savo ribose nupaisomi
    if (testid>1) {
      for (i in 2:c(testid)){
        print(c("brai?om r?m?", i-1))
        lines(x=rims[[i]][,1],y=rims[[i]][,2],col=1,lwd=2)
      }}
    #testavimui pjuviai paruosiami
    t<-per.taskai(margins,n=divisions)
    points(t[[3]],pch=19,col="purple")
    linijos<-line2(t[[1]],t[[2]],polygon = margins)
    maxdif<- original.qual
    print(maxdif)
    any.split<-numeric()
    maxid<-0
    if (nrow(linijos)!=0){
      points(linijos[,c(1:2)],col="yellow",pch=19)
      points(linijos[,c(3:4)],col="red",pch=19)
      #pjaustymo ir testavimo ciklas
      {
        for (i in 1:nrow(linijos)){
          print('testuojamas padalinimas Nr.:')
          print(i)
          virs.h<-filter.reorg.poly(polygon = t[[3]][-nrow(t[[3]]),],poli.side = T, min.x.id=linijos[i,6],max.x.id=linijos[i,7],b=linijos[i,5])
          virs<-close.poly(virs.h)
          po.h<-filter.reorg.poly(polygon = t[[3]][-nrow(t[[3]]),],poli.side = F, min.x.id=linijos[i,6],max.x.id=linijos[i,7],b=linijos[i,5])
          po<-close.poly(po.h)

          # padalinami duomenys i dvi dalis pagal pjuvio koordinates
          Puses<-list(get.data(po,samp.dat),get.data(virs,samp.dat))
          if (all(c(length(Puses[[1]]$rusis),length(Puses[[2]]$rusis))>N.cond)){
            SpjuvioI<-abs(polyarea(x=virs[,1],y=virs[,2]))
            SpjuvioII<-abs(polyarea(x=po[,1],y=po[,2]))
            if (SpjuvioI>S.cond&SpjuvioII>S.cond){
              #nupiesiam padalinima ir paskaiciuojam kokybe
              Skirtumas<-alfa(Puses[[1]]$rusis,Puses[[2]]$rusis)
              any.split<-c(any.split,Skirtumas)
              #Paskaiciuojam plotus padalintu bloku
              if (Skirtumas > maxdif){
                #Jei padalinimas patenkina minimalias saligas ir yra geresnis nei pries tai - pasizymim ir issisaugom ji
                maxdif<-Skirtumas
                maxid<-i
                print(c('max Skirtumas=',Skirtumas))
              }
            }
          }
        }
      }
      if(length(any.split)>0){
        performance<-mean(any.split)
        #jei performance 0 ir visu padalinimu performance 0, tai padalinimo nera - idealiai atskirtas plotas
        if(all(any.split==0)){
          maxid<-0
        }
      } else {
        performance<-original.qual # negalejom ivertinti ne vieno padalinimo, taigi performance lygu max.
      }

      blokai<<-rbind(blokai,data.frame(performance=1-performance/original.qual,
                                      iteration=iteration,saknis=root))
      print(c("blokai: ", blokai))
      print(c("performance: ", performance))
    }
    # duomenu saugojimas
    #Jei rastas tinkamas padalinimas - ieskom geriausios padalinimo kreives,
    #issaugom duomenis ir ziurim ar galima skaidyti toliau
    if (maxid>0){
      best.curve<-curvial.split(poly.x=t[[3]]$x.polio,poly.y=t[[3]]$y.polio,
                                min.x.id = linijos[maxid,6],max.x.id = linijos[maxid,7],b=linijos[maxid,5],samp.dat,knot.density.X=knot.density.X,
                                knot.density.Y=knot.density.Y,N.cond,S.cond,iteracija=curve.iterations,correction.term=correction.term,original.qual)
      lines(best.curve[[1]],col=2,lwd=3)
      if ((1-(max(best.curve[[2]],maxdif)/original.qual))<C.cond ){
        maxid<-0}}
    if (maxid>0){
      if(best.curve[[2]]>maxdif) {
        splits<-do.call(c,list(splits,list(data.frame(x=best.curve[[1]]$x,y=best.curve[[1]]$y))))
        split.performance<-c(split.performance,best.curve[[2]])
        split.reliability<-c(split.reliability,(best.curve[[2]]-performance)/sd(any.split))
        #dalinam duomenis padalinimo kreive
        #reikia sukurti poligonus du ir nufiltruoti duomenis  - galima padaryti geriau
        virsus.h<-filter.reorg.poly(polygon = data.frame(xp=t[[3]][,1][-nrow(t[[3]])],yp=t[[3]][,2][-nrow(t[[3]])]),min.x.id =linijos[maxid,6],
                                    max.x.id =linijos[maxid,7],poli.side = T,b= linijos[maxid,5])
        O.poli<-close.poly(split.poly = virsus.h,split.line.x = best.curve[[1]]$x, split.line.y = best.curve[[1]]$y)
        O<-get.data(O.poli,samp.dat)
        lines(O.poli,col=5,lwd=4)

        apacia.h<-filter.reorg.poly(polygon = data.frame(xp=t[[3]][,1][-nrow(t[[3]])],yp=t[[3]][,2][-nrow(t[[3]])]),min.x.id =linijos[maxid,6],max.x.id =linijos[maxid,7],poli.side = F,b=linijos[maxid,5])
        OO.poli<-close.poly(split.poly = apacia.h,split.line.x = best.curve[[1]]$x, split.line.y = best.curve[[1]]$y)
        lines(OO.poli,col=5,lwd=4)

        OO<-get.data(OO.poli,samp.dat)

        #issaugom duomenis padalinimo
        ribs<-list(O.poli,OO.poli)
      } else {
        #issaugom duomenis padalinimo
        splits<-do.call(c,list(splits,list(data.frame(x=as.numeric(c(linijos[maxid,c(1,3)])),y=as.numeric(c(linijos[maxid,c(2,4)]))))))
        split.performance<-c(split.performance,maxdif)
        split.reliability<-c(split.reliability,(maxdif-performance)/sd(any.split))
        #padalinam duomenis pagal pjuvi
        virs<-close.poly(split.line.x = as.numeric(linijos[maxid,c(1,3)]),split.line.y = as.numeric(linijos[maxid,c(2,4)]),
                         split.poly =filter.reorg.poly(polygon = data.frame(xp=t[[3]][,1][-nrow(t[[3]])],yp=t[[3]][,2][-nrow(t[[3]])]),
                                                       min.x.id =linijos[maxid,6],max.x.id =linijos[maxid,7],poli.side = T,b=linijos[maxid,5]))

        po<-close.poly(split.line.x = as.numeric(linijos[maxid,c(1,3)]),split.line.y = as.numeric(linijos[maxid,c(2,4)]),
                       split.poly =filter.reorg.poly(polygon = data.frame(xp=t[[3]][,1][-nrow(t[[3]])],yp=t[[3]][,2][-nrow(t[[3]])]),
                                                     min.x.id =linijos[maxid,6],max.x.id =linijos[maxid,7],poli.side = F,b=linijos[maxid,5]))

        O<-get.data(virs,samp.dat)
        OO<-get.data(po,samp.dat)

        #paryskinam atrinkta pjuvi
        print(c('total max Skirtumas =', maxdif))

        lines(linijos[maxid,c(1,3)],linijos[maxid,c(2,4)],col=2)

        # sukuriam koordinates ribines zemu koordinaciu plotui ir
        # bandom skaidyti ji. Jei pavyksta suskaidyti issaugom updatintus rezus, jei ne istrinam gogolis
        # ir ribines koordinates neegzistuojancio pjuvio
        ribs<-list(virs,po)
      }
      #maisom duomenis ir vertinam aptiktu erdviniu strukturu patimuma
      if (null.models){
        set.seed(seed.t)
        test.samp.dat<-sim.testdat(samp.dat,test.n)
        pseudo.kokybe<-numeric(test.n)
        for (a in 1:test.n){
          I.dat<-get.data(ribs[[1]],test.samp.dat[[a]])
          II.dat<-get.data(ribs[[2]],test.samp.dat[[a]])
          pseudo.kokybe[a]<-alfa(I.dat[,1],II.dat[,1])
        }
        checks<-do.call(c,list(checks,list(pseudo.kokybe)))
        split.reliability2<-c(split.reliability2,sum(last(split.performance)<pseudo.kokybe)/test.n)
      }
      n.splits<-c(n.splits,length(any.split))
      ave.split.abE<-c(ave.split.abE,performance)
      original.quality<<-c(original.quality,original.qual)
      #
      rims<-do.call(c,list(rims,ribs[1]))
      lines(ribs[[1]],col="purple")
      gogolis<-spatial_div(O,
                           rims = rims,divisions = divisions,
                           method=method,variation=variation,metric=metric,criteria=criteria,
                           C.cond=C.cond,E.cond=E.cond,N.cond=N.cond,S.cond=S.cond,
                           root = iteration,knot.density.X = knot.density.X,
                           knot.density.Y = knot.density.Y,curve.iterations = curve.iterations,correction.term=correction.term,
                           split.reliability = split.reliability,splits=splits,
                           split.performance = split.performance,split.reliability2=split.reliability2,checks=checks,
                           null.models=null.models,n.splits = n.splits,seed.t=seed.t,ave.split.abE = ave.split.abE,test.n = test.n)
      print(paste("griztam i", testid, "padalinima [po mazu koord bloko]", sep=" "))

      if (length(gogolis)==2){

        blokai<-gogolis[[1]]
        rims<-gogolis[[2]]
        rm(gogolis)
      }

      # Skaidom antra bloka
      #jei egzistuoja gogolis (updatinti duomenys) tuomet sukuriam ribines koordinates ir lipdom prie
      #gogolis masyvo. Duotu koordinaciu ribose ir bandom ieskoti pjuvio bei toliau updatinti duomenis.
      #Jei pavyksta rasti pjuvi, updatinam duomenis, jei ne trinam null bobolis ir priklijuotas koo-
      #rdinates.
      if (exists("gogolis")){
        gogolis[[2]]<-do.call(c,list(gogolis[[2]],ribs[2]))
        lines(ribs[[2]],col=2)
        bobolis<-spatial_div(OO,divisions = divisions,
                             rims = gogolis[[2]],
                             method=method,variation=variation,metric=metric,criteria=criteria,
                             C.cond=C.cond,E.cond=E.cond,N.cond=N.cond,S.cond=S.cond,
                             root=iteration,
                             knot.density.X = knot.density.X,knot.density.Y = knot.density.Y,curve.iterations = curve.iterations,
                             correction.term=correction.term,
                             split.reliability = gogolis[[4]],null.models = null.models,
                             splits = gogolis[[1]],split.performance=gogolis[[5]],split.reliability2=gogolis[[7]],
                             checks=gogolis[[8]],n.splits = gogolis[[9]],seed.t=seed.t,ave.split.abE = gogolis[[10]],test.n = test.n)
        print(paste("griztam i", testid, "padalinima [po aukstu koord bloko (gogolis exists)]", sep=" "))

        if (length(bobolis)==2){

          gogolis[[3]]<-bobolis[[1]]
          gogolis[[2]]<-bobolis[[2]]
          rm(bobolis)
        }

      } else{
        #jei pirmo bloko skaidymas nebuvo sekmingas, sukuriam ribines koordinates naujam pjuviui
        # ir bandom skaidyti.
        rims<-do.call(c,list(rims,ribs[2]))
        lines(ribs[[2]],col=2)

        bobolis<-spatial_div(OO,
                             rims = rims,divisions = divisions,
                             method=method,variation=variation,metric=metric,criteria=criteria,
                             C.cond=C.cond,E.cond=E.cond,N.cond=N.cond,S.cond=S.cond,
                             root = iteration,knot.density.X = knot.density.X,
                             knot.density.Y = knot.density.Y,curve.iterations = curve.iterations,correction.term=correction.term,
                             split.reliability = split.reliability,splits=
                               splits,split.performance=split.performance,split.reliability2=split.reliability2,checks=checks,
                             null.models = null.models,n.splits = n.splits,seed.t=seed.t,ave.split.abE = ave.split.abE,test.n = test.n)
        print(paste("griztam i", testid, "padalinima [po aukstu koord bloko (no gogolis)]", sep=" "))
        # Jei bobolis updatino duomenis ji paliekam, jei ne istrinam, kartu pasalindami koordinates

        if (length(bobolis)==2){
          blokai<-bobolis[[1]]
          rims<-bobolis[[2]]
          rm(bobolis)
        }

      }
      #jei duomenys buvo sukurti - siunciami nauji, jei ne - siunciami su naujo padalinimo info,
      # be koordinaciu nauju plotu skaidymui
      if (exists("bobolis")){
        print("grazinami updatinti duomenis po antro skaidymo su naujom ribom [?+-pirmas]")

        return(bobolis)} else{
          if (exists("gogolis")) {
            print("grazinami updatinti duomenis po pirmo skaidymo su naujom ribom, antras nesekmingas")

            return(gogolis)} else{
              print("pasiekta nebesiskaidymo riba, grazinami updatinti duomenis, be nauju ribiniu")
              return(list(splits,rims,1,split.reliability,split.performance,1,split.reliability2,checks,n.splits,ave.split.abE))
            }
        }
    } else{
      #Jei tinkamo padalino nerasta, grizta tuscias masyvas
      if (testid>1){
        print("tinkamo padalinimo nerasta, grizta tuscias masyvas NULL")
        return(list(1,rims))
      } else{
        print("nebuvo skaldomu bloku")
        return(list(1,rims))}
    }
  }
  environment(spatial_div) <- environment()

  rezas<-spatial_div(data,method=method,variation=variation,metric=metric,criteria=criteria,
                     C.cond=C.cond,E.cond=E.cond,N.cond=N.cond,S.cond=S.cond,rims = rims,divisions = divisions,
                     root=2,knot.density.X=knot.density.X,
                     knot.density.Y=knot.density.Y,curve.iterations=curve.iterations,correction.term=correction.term,
                     split.reliability = split.reliability,splits=splits,split.performance = split.performance,
                     split.reliability2=split.reliability2,checks=checks,null.models = null.models,n.splits = n.splits,seed.t=seed.t,
                     ave.split.abE = ave.split.abE,test.n = test.n)


  if (null.models==F){
    rezas[[7]]<-rep(NaN,length(rezas[[4]]))
  }

  Signif <- symnum(rezas[[7]], corr = FALSE, na = FALSE,
                   cutpoints = c(0 ,0.001,0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))


  rezas <- structure(list(
    split.lines = rezas[[1]],
    boundaries = rezas[[2]],
    block.stats = blokai[-1,],
    split.stats = data.frame(
      n.splits = rezas[[9]],
      z.score = round(rezas[[4]],2),
      ave.abE = round(-rezas[[10]],2),
      split.abE = round(-rezas[[5]],2),
      parent.2E = round(-original.quality,2),
      delta.E = round(rezas[[5]] -  original.quality,2),
      p_value = rezas[[7]],
      signif. = format(Signif)
    ),
    null.m.st= rezas[[8]]
  ),
  class = "spdiv"
  )
  print.spdiv <- function(x){
    cat("\n","Information about splits:", "\n","\n")
    print(x[[4]])
    cat("\n", "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    invisible(x)
  }
  return(print.spdiv(rezas))
}

o<-spatial.analysis(data = data,C.cond = 0.05,N.cond = 1000*n/n.si,S.cond = 0.1,E.cond = 0,divisions = 15,knot.density.X = 6,
                    knot.density.Y = 30,curve.iterations = 3,correction.term=0.1,null.models = F,seed.t = 2,test.n = 100)
plot(o,data)
o$block.stats
performance	iteracija	saknis
2	0.152951876163506	2	2
3	0	                3	2
4	0.175787975012891	4	2
5	0.172868580412207	5	4
6	0	               	6	5
7	0	               	7	5
8	0	               	8	4



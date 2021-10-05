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
  spatial_div<-function(samp.dat, root=2){
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
    perim_pts<-.perimeter_pts(polygon = margins,n.pts = divisions)
    points(perim_pts[[2]],pch=19,col="purple")
    pairs_pts<-.pair_pts(perim_pts[[1]],polygon = margins)
    maxdif<- original.qual
    print(maxdif)
    any.split<-numeric()
    maxid<-0
    if (nrow(pairs_pts)!=0){
      points(pairs_pts[,c(1:2)],col="yellow",pch=19)
      points(pairs_pts[,c(3:4)],col="red",pch=19)
      #pjaustymo ir testavimo ciklas
      {
        for (i in 1:nrow(pairs_pts)){
          print('testuojamas padalinimas Nr.:')
          print(i)

          virs <- .close_poly(.split_poly(polygon = perim_pts[[2]], min_id = 1,
                              split_ids = pairs_pts[i,c(6:7)],
                              trivial_side = TRUE,poli_side = TRUE))
          po <- .close_poly(.split_poly(polygon = perim_pts[[2]], min_id = 1,
                            split_ids = pairs_pts[i,c(6:7)],
                            trivial_side = TRUE,poli_side = FALSE))

          # padalinami duomenys i dvi dalis pagal pjuvio koordinates
          Puses <- list(.get_data(po,samp.dat),.get_data(virs,samp.dat))
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
      print(perim_pts)
      print(pairs_pts)
      best.curve<-curvial.split(poly.x=perim_pts[[2]]$x.poly,poly.y=perim_pts[[2]]$y.poly,
                                min.x.id = pairs_pts[maxid,6],max.x.id = pairs_pts[maxid,7],b=pairs_pts[maxid,5],samp.dat,knot.density.X=knot.density.X,
                                knot.density.Y=knot.density.Y,N.cond,S.cond,iteracija=curve.iterations,correction.term=correction.term,original.qual)
      lines(best.curve[[1]],col=2,lwd=3)
      if ((1-(max(best.curve[[2]],maxdif)/original.qual))<C.cond ){
        maxid<-0}}
    if (maxid>0){
      if(best.curve[[2]]>maxdif) {
        splits<<-do.call(c,list(splits,list(data.frame(x=best.curve[[1]]$x,y=best.curve[[1]]$y))))
        split.performance<<-do.call(c,list(split.performance,best.curve[[2]]))
        split.reliability<<-do.call(c,list(split.reliability,(best.curve[[2]]-performance)/sd(any.split)))
        #dalinam duomenis padalinimo kreive
        #reikia sukurti poligonus du ir nufiltruoti duomenis  - galima padaryti geriau
        virsus.h<-filter.reorg.poly(polygon = data.frame(xp=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],yp=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),min.x.id =pairs_pts[maxid,6],
                                    max.x.id =pairs_pts[maxid,7],poli.side = T,b= pairs_pts[maxid,5])
        O.poli<-close.poly(split.poly = virsus.h,split.line.x = best.curve[[1]]$x, split.line.y = best.curve[[1]]$y)
        O<-get.data(O.poli,samp.dat)
        lines(O.poli,col=5,lwd=4)

        apacia.h<-filter.reorg.poly(polygon = data.frame(xp=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],yp=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),min.x.id =pairs_pts[maxid,6],max.x.id =pairs_pts[maxid,7],poli.side = F,b=pairs_pts[maxid,5])
        OO.poli<-close.poly(split.poly = apacia.h,split.line.x = best.curve[[1]]$x, split.line.y = best.curve[[1]]$y)
        lines(OO.poli,col=5,lwd=4)

        OO<-get.data(OO.poli,samp.dat)

        #issaugom duomenis padalinimo
        ribs<-list(O.poli,OO.poli)
      } else {
        #issaugom duomenis padalinimo
        splits<<-do.call(c,list(splits,list(data.frame(x=as.numeric(c(pairs_pts[maxid,c(1,3)])),y=as.numeric(c(pairs_pts[maxid,c(2,4)]))))))
        split.performance<<-do.call(c,list(split.performance,maxdif))
        split.reliability<<-do.call(c,list(split.reliability,(maxdif-performance)/sd(any.split)))
        #padalinam duomenis pagal pjuvi
        virs<-close.poly(split.line.x = as.numeric(pairs_pts[maxid,c(1,3)]),split.line.y = as.numeric(pairs_pts[maxid,c(2,4)]),
                         split.poly =filter.reorg.poly(polygon = data.frame(xp=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],yp=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
                                                       min.x.id =pairs_pts[maxid,6],max.x.id =pairs_pts[maxid,7],poli.side = T,b=pairs_pts[maxid,5]))

        po<-close.poly(split.line.x = as.numeric(pairs_pts[maxid,c(1,3)]),split.line.y = as.numeric(pairs_pts[maxid,c(2,4)]),
                       split.poly =filter.reorg.poly(polygon = data.frame(xp=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],yp=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
                                                     min.x.id =pairs_pts[maxid,6],max.x.id =pairs_pts[maxid,7],poli.side = F,b=pairs_pts[maxid,5]))

        O<-get.data(virs,samp.dat)
        OO<-get.data(po,samp.dat)

        #paryskinam atrinkta pjuvi
        print(c('total max Skirtumas =', maxdif))

        lines(pairs_pts[maxid,c(1,3)],pairs_pts[maxid,c(2,4)],col=2)

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
        checks<<-do.call(c,list(checks,list(pseudo.kokybe)))
        split.reliability2<<-do.call(c,list(split.reliability2,sum(last(split.performance)<pseudo.kokybe)/test.n))
      }
      n.splits<<-do.call(c,list(n.splits,length(any.split)))
      ave.split.abE<<-do.call(c,list(ave.split.abE,performance))
      original.quality<<-do.call(c,list(original.quality,original.qual))
      #
      rims<<-do.call(c,list(rims,ribs[1]))
      lines(ribs[[1]],col="purple")
      spatial_div(O,root = iteration)
      print(paste("griztam i", testid, "padalinima [po mazu koord bloko]", sep=" "))


      # Skaidom antra bloka
      #jei egzistuoja gogolis (updatinti duomenys) tuomet sukuriam ribines koordinates ir lipdom prie
      #gogolis masyvo. Duotu koordinaciu ribose ir bandom ieskoti pjuvio bei toliau updatinti duomenis.
      #Jei pavyksta rasti pjuvi, updatinam duomenis, jei ne trinam null bobolis ir priklijuotas koo-
      #rdinates.

        rims<<-do.call(c,list(rims,ribs[2]))
        lines(ribs[[2]],col=2)
        spatial_div(OO,root=iteration)
        print(paste("griztam i", testid, "padalinima [po aukstu koord bloko (gogolis exists)]", sep=" "))

    } else{
      #Jei tinkamo padalino nerasta, grizta tuscias masyvas
      if (testid>1){
        print("tinkamo padalinimo nerasta, grizta tuscias masyvas NULL")
      } else{
        print("nebuvo skaldomu bloku")
}
    }
  }
  environment(spatial_div) <- environment()

  spatial_div(data,root=2)


  if (null.models==F){
    split.reliability2<-rep(NaN,length(n.splits))
  }

  Signif <- symnum(split.reliability2, corr = FALSE, na = FALSE,
                   cutpoints = c(0 ,0.001,0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))


  rezas <- structure(list(
    split.lines = splits,
    boundaries = rims,
    block.stats = blokai[-1,],
    split.stats = data.frame(
      n.splits = n.splits,
      z.score = round(split.reliability,2),
      ave.abE = round(-ave.split.abE,2),
      split.abE = round(-split.performance,2),
      parent.2E = round(-original.quality,2),
      delta.E = round(split.performance -  original.quality,2),
      p_value = split.reliability2,
      signif. = format(Signif)
    ),
    null.m.st= checks
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



#' Main recursive hespdiv inner helper function
#'
#' @description  During each recursive iteration this function fits straight and
#' curvi-linear split-lines, then it uses the better of the two to separate
#' data in space, saves information about the polygon and split-lines,
#' and calls itself again using the extracted data samples.
#' @param samp.dat a spatial subset of \code{data} that lies entirely within
#' some polygon, produced using split-line subdivisions of original study area.
#' At first iteration, while no subdivisions are established,
#' \code{samp.dat = data}.
#' @param samp.xy a spatial subset of \code{xy.dat}
#' @param root.id An id of recursive iteration that produced the \code{samp.dat}.
#' An id of parent-polygon.
#' @return No return. Function updates variables in \code{hespdiv} environment.
#' @author Liudas Daumantas
#' @note Function inherits the environment of \code{hespdiv} function.
#' @importFrom pracma poly_center
#' @noRd

.spatial_div <- function(samp.dat, samp.xy, root.id){

  #testuojamas plotas
  testid <- length(rims)

  assign(x = "poly.obj" ,
         value = do.call(c,list(poly.obj,list(generalize.f(samp.dat)))),
         envir = e)
  perim_pts <- .perimeter_pts(polygon = rims[[testid]],n.pts = n.split.pts)
  #if (testid == 28){
  #trace.object <- "both"
  # trace.level <- "all"
  #}
  .visualise_splits.start(what = trace.object,
                          pnts.col, xy.dat, rims,perim_pts,
                          testid)

  #testavimui pjuviai paruosiami

  environment(.dif_fun) <- e

  pairs_pts <- .pair_pts(perim_pts[[1]],polygon = rims[[testid]])
  maxdif <- lower.Q.crit # first split minimum quality. P.crit
  any.split <- numeric()
  maxid <- 0

  if (nrow(pairs_pts) != 0){
    #pjaustymo ir testavimo ciklas
    for (i in 1:nrow(pairs_pts)){
      virs <- .close_poly(
        open.poly =
          .split_poly(
            polygon = perim_pts[[2]],
            min_id = 1,
            split_ids = as.numeric(pairs_pts[i,c(6:7)]),
            trivial_side = TRUE,
            poli_side = TRUE
          ))
      po <- .close_poly(
        open.poly =
          .split_poly(
            polygon = perim_pts[[2]],
            min_id = 1,
            split_ids = as.numeric(pairs_pts[i,6:7]),
            trivial_side = TRUE,
            poli_side = FALSE
          ))

      # padalinami duomenys i dvi dalis pagal pjuvio koordinates
      id1 <- .get_ids(po, samp.xy,first.p, data.frame(
          x = unlist(pairs_pts[i,c(1,3)]), y = unlist(pairs_pts[i,c(2,4)])) )
      id2 <- .get_ids(virs, samp.xy,first.p, data.frame(
          x = unlist(pairs_pts[i,c(1,3)]), y = unlist(pairs_pts[i,c(2,4)])) )


      .visualise_splits.try_straight(what = trace.object, level = trace.level,
                                     pairs_pts, i)

      if (all(c(length(id1), length(id2)) > N.crit)){
        if (S.crit > 0){
          SpjuvioI <- abs(pracma::polyarea(x=virs[,1],y=virs[,2]))
          SpjuvioII <- abs(pracma::polyarea(x=po[,1],y=po[,2]))
          cond <- SpjuvioI > S.cond & SpjuvioII > S.cond
        } else {
          cond <- TRUE
        }
        if (cond) {
          #nupiesiam padalinima ir paskaiciuojam kokybe
          environment(.dif_fun) <- environment()
          Skirtumas <- .dif_fun(.slicer(samp.dat,id1), .slicer(samp.dat,id2))
          any.split <- c(any.split,Skirtumas)
          #Paskaiciuojam plotus padalintu bloku
          if (.comp(Skirtumas,maxdif) ){

            .visualise_splits.good_straight(what = trace.object,
                                            level = trace.level,
                                            pairs_pts, Skirtumas, i,
                                            maxid)

            #Jei padalinimas patenkina minimalias saligas ir yra
            # geresnis nei pries tai - pasizymim ir issisaugom ji
            maxdif <- Skirtumas
            maxid <- i
          } else {
            message <- paste0("Poor split quality.\n","Obtained: ",
                              round(Skirtumas,2),
                              "\nRequired: ",c.sign, round(maxdif,2))
          }
        } else {
          message <- paste0("One of the areas was too small.\n","Obtained: ",
                            round(SpjuvioI,2), ' and ', round(SpjuvioII,2),
                            "\nRequired: >",S.cond)
        }
      } else {
        message <- paste0("Not enough observations in one of the areas.",
                          "\nObtained: ", length(id1), ' and ', length(id2),
                          "\nRequired: >",N.crit)
      }
      if (maxid != i) {
        .visualise_splits.bad_straight(what = trace.object,
                                       level = trace.level,
                                       pairs_pts, message, i)
      }
    }
    assign(x = "str.split.quals" ,
           value = do.call(c,list(str.split.quals,list(any.split))),
           envir = e)

    if(length(any.split) > 0){

      mean.dif <- mean(any.split) # mean spatial heterogeneity
      sd.dif <- sd(any.split) # NA if 1 split
      # sd(any.split) -> anisotropy of heterogeneity
      # if all are equal, then complete randomness: no anisotropy - no splits present
      if( all(any.split==any.split[1]) & length(any.split) > 1  ){
        warning(paste(c(
          "All tested splits were of equal quality in "
          , testid, " iteration. A random split were selected as it meets",
          " all provided criteria")))
      }
    } else {
      mean.dif <- NA # negalejom ivertinti ne vieno padalinimo, taigi performance lygu max.
      sd.dif <- NA # P.crit
    }
    assign(x = "poly.info" ,value = rbind(poly.info, data.frame(
      root.id = root.id,
      n.splits = length(any.split),
      n.obs = nrow(samp.xy),
      mean = mean.dif, # NA if 0
      sd = sd.dif, # NA if 1 or 0 split
      str.best = maxdif,
      str.z.score = (maxdif - mean.dif) / sd.dif, # NA if 1 or 0 split
      has.split = FALSE,
      is.curve = FALSE,
      crv.best = NA,
      crv.z.score = NA,
      c.improv = NA)
    )
    , envir = e)
    # duomenu saugojimas
    #Jei rastas tinkamas padalinimas - ieskom geriausios padalinimo kreives,
    #issaugom duomenis ir ziurim ar galima skaidyti toliau
    if (maxid>0) {
      if ( c.splits) {
        .visualise_splits.best_straight(what = trace.object,
                                        pairs_pts, maxid, maxdif)
        environment(.curvial_split) <- environment()
        best.curve <- .curvial_split(

          poly.x = perim_pts[[2]]$x.poly,
          poly.y = perim_pts[[2]]$y.poly,
          min.x.id = pairs_pts[maxid,6],
          max.x.id = pairs_pts[maxid,7],
          b = pairs_pts[maxid,5],
          samp.dat = samp.dat,
          c.X.knots = c.X.knots,
          c.Y.knots = c.Y.knots,
          N.cond = N.crit,
          S.cond = S.cond,
          c.max.iter.no = c.max.iter.no,
          c.fast.optim = c.fast.optim,
          c.corr.term = c.corr.term,
          trace.object = trace.object,
          trace.level = trace.level,
          pnts.col = pnts.col,
          straight.qual = maxdif

        )
        poly.info[testid,"crv.best"] <- best.curve[[2]]
        poly.info[testid,"crv.z.score"] <- (best.curve[[2]] - mean.dif) / sd.dif
        poly.info[testid,"c.improv"] <- abs(best.curve[[2]] - maxdif)
        if(.comp(best.curve[[2]],maxdif) &
           poly.info[testid,"c.improv"] >= C.crit.improv){
          poly.info[testid,"is.curve"] <- TRUE
          maxdif <- best.curve[[2]]
        }
        if (!.comp(maxdif, upper.Q.crit) ){
          maxid <- 0
        }
      } else {
        if (!.comp(maxdif,upper.Q.crit) ){
          maxid <- 0
        }
      }
    }

    if (maxid > 0){ # save the split and perform new splits if TRUE

      poly.info[testid,"has.split"] <- TRUE
      assign(x = "poly.info" ,value = poly.info, envir = e)

      if(poly.info[testid,"is.curve"]) {
        best.splitl <- data.frame(x = best.curve[[1]]$x, y = best.curve[[1]]$y)
      } else {
        best.splitl <- data.frame(x=as.numeric(c(pairs_pts[maxid,c(1,3)])),
                                  y=as.numeric(c(pairs_pts[maxid,c(2,4)])))
      }
      .visualise_splits.best_split(what = trace.object,
                                   best.splitl, maxdif)

      assign(x = "splits" ,value = do.call(c,list(splits,list(best.splitl))),
             envir = e)
      # z-score of performance
      #dalinam duomenis padalinimo kreive
      #reikia sukurti poligonus du ir nufiltruoti duomenis  - galima padaryti geriau
      up.pol <- .close_poly(
        open.poly =
          .split_poly(
            polygon = data.frame(x=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],
                                 y=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
            split_ids = as.numeric(pairs_pts[maxid,6:7]),
            min_id = 1,
            trivial_side = TRUE,
            poli_side = TRUE
          ),
        close.line = best.splitl
      )
      up.ids <- .get_ids(up.pol,samp.xy,first.p, data.frame(
        x = unlist(pairs_pts[maxid,c(1,3)]), y = unlist(pairs_pts[maxid,c(2,4)])))
      up.xy <- .slicer.table(samp.xy,up.ids)
      up.dat <- .slicer(samp.dat,up.ids)

      do.pol <- .close_poly(
        open.poly =
          .split_poly(
            polygon = data.frame(x=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],
                                 y=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
            split_ids = as.numeric(pairs_pts[maxid,6:7]),
            min_id = 1,
            trivial_side = TRUE,
            poli_side = FALSE
          ),
        close.line = best.splitl
      )
      do.ids <- .get_ids(do.pol,samp.xy,first.p, data.frame(
        x = unlist(pairs_pts[maxid,c(1,3)]), y = unlist(pairs_pts[maxid,c(2,4)])))
      do.xy <- .slicer.table(samp.xy,do.ids)
      do.dat <- .slicer(samp.dat,do.ids)

      #issaugom duomenis padalinimo
      # ribs <- list(up.pol,do.pol)

      #maisom duomenis ir vertinam aptiktu erdviniu strukturu patimuma
      if (n.m.test){

        set.seed(n.m.seed)
        test1 <- .sp.n.m(samp.dat,ribs,n.m.N,n.m.keep,type = 1)
        test2 <- .sp.n.m(samp.dat,ribs,n.m.N,n.m.keep,type = 2)

        if (n.m.keep){
          assign(x = "n.m.sims1",value = do.call(c,list(n.m.sims1,list(test1[[2]]))),
                 envir = e)
          assign(x = "n.m.sims2",value = do.call(c,list(n.m.sims2,list(test2[[2]]))),
                 envir = e)
        }
        assign(x = "sim2.difs",value = do.call(c,list(sim2.difs,list(test2[[1]]))),
               envir = e)
        assign(x = "sim1.difs",value = do.call(c,list(sim1.difs,list(test1[[1]]))),
               envir = e)

        assign(x = "p.val1", value =
                 do.call(c,list(
                   p.val1,
                   sum(maxdif < test1[[1]])/n.m.N
                 )),
               envir = e) # kvantilis
        # kvantilio reiksme - empirine p verte
        # jei gausinis skirstinys,tai:
        # 1 - qnorm(sum(last(split.quality),mean(pseudo.kokybe),sd(pseudo.kokybe)) duotu teorine p-verte
        assign(x = "p.val2", value =
                 do.call(c,list(
                   p.val,
                   sum(maxdif<test2[[1]])/n.m.N
                 )),
               envir = e)
      }

      assign(x = "rims" ,value = do.call(c,list(rims,list(up.pol))) ,envir = e)

      .spatial_div(up.dat,up.xy, root.id = testid)

      print(paste("griztam i", testid, "padalinima [po mazu koord bloko]", sep=" "))


      # Skaidom antra bloka
      #jei egzistuoja gogolis (updatinti duomenys) tuomet sukuriam ribines koordinates ir lipdom prie
      #gogolis masyvo. Duotu koordinaciu ribose ir bandom ieskoti pjuvio bei toliau updatinti duomenis.
      #Jei pavyksta rasti pjuvi, updatinam duomenis, jei ne trinam null bobolis ir priklijuotas koo-
      #rdinates.
      assign(x = "rims" ,value = do.call(c,list(rims,list(do.pol))) ,envir = e)
      .spatial_div(do.dat,do.xy, root.id = testid)
      print(paste("griztam i", testid, "padalinima [po aukstu koord bloko (gogolis exists)]", sep=" "))

    }
    } else{
      #Jei tinkamo padalino nerasta, grizta tuscias masyvas
      if (testid>1){
        if (!is.null(trace.object))
          print("There were no suitable points on the polygon perimeter
            to generate straight-split lines")
      } else{
        if (!is.null(trace.object))
          print("There were no suitable points on the polygon perimeter
            to generate straight-split lines")
      }
    }
}



#' Calculate difference between two data sets
#'
#' @description  Function combines generalize.f and compare.f functions, to
#' estimate the difference between two data sets obtained from opposite sides
#' of a split-line. This difference reflects the ability of a split-line to
#' spatially separate data. Thus, it is considered to reflect the quality of a
#' split-line.
#' @param dat1 a spatial subset of \code{samp.dat} that is located at opposite
#' side of a split-line than \code{dat2} is.
#' \code{samp.dat = data}.
#' @param dat2 a spatial subset of \code{samp.dat} that is located at opposite
#' side of a split-line than \code{dat1} is.
#' @return Numeric value, expressing the difference between two data sets and
#' the quality of a split-line.
#' @author Liudas Daumantas
#' @note Function inherits the environment of \code{hespdiv} function. Also,
#' it forces this environment to be inherited by generalize.f and compare.f
#' functions.
#' @noRd
.dif_fun <- function(dat1,dat2) {
  environment(generalize.f) <- environment()
  environment(compare.f) <- environment()
  compare.f( generalize.f(dat1), generalize.f(dat2) )
}




#' Test the split-line with null model simulations
#'
#' @description  Function simulates spatial null models of the provided data and
#' then checks the performance of the established split-line fitted to the
#' simulated data.
#' @param data a samp.dat from .spatial_div function.
#' @param ribs a list of two data frames of polygons, established by the
#' fitted, "best" split-line.
#' @param n integer - number of data simulations to perform.
#' @param n.m.keep logical - should the produced simulations be kept
#' @param type type of null model simulations:
#' 1 - completely random: toroidal rotation of individual data points (the only
#' things maintained is polygon shape and number of points, and points
#' themselves),
#' 2 - random perspective of data points: toroidal rotation of all data points
#' (micro and macro spatial inter-relations of points are maintained)
#' Other possible n.m. types:
#' 3 - random macro inter-relations of data points: toroidal rotation of
#' data point groups, that lie in the same quadrant of polygon extent (micro
#' inter-relations of points are maintained) # quadrant number would be required
#' as an argument.
#' 4 - random micro inter-relations of data points: jitter of data
#' points. The amount of jittery equal to points' distance to the xth
#' nearest neighbor (density surface (intensity, first moment) as well as macro
#' inter-relations of points are maintained) # xth nearest neighbour would be
#' required as an argument.
#' 5?? - random micro and macro inter-relations of data points: rotation of
#' quadrants and jitter of the points (maintained are point inter-relations that
#' are occur between the given micro and macro scales.)
#' @return list of one or two elements: 1 - vector of estimated performances of
#' the split-line for each of null model simulation. 2 -  list of data frames of
#' simulated data (null models), returned only when n.m.keep is TRUE.
#' @author Liudas Daumantas
#' @noRd
.sp.n.m <- function(data,ribs,n,n.m.keep,type){
  if (type == 1) {
    N <- nrow(data)
  } else {
    N <- 1
  }
  mirror.data <- data
  pseudo.kokybe <- numeric(n.m.N)
  if (n.m.keep){
    sim <- list()
  }
  for (i in 1:n){
    x.shift <- runif(N,min=0,max=dist(range(data$x)))
    y.shift <- runif(N,min=0,max=dist(range(data$y)))
    testx <- data$x+x.shift
    testy <- data$y+y.shift
    tx <- case_when(testx>max(data$x) ~ testx-max(data$x)+min(data$x),
                    TRUE ~ testx)
    ty <- case_when(testy>max(data$y) ~ testy-max(data$y)+min(data$y),
                    TRUE ~ testy)

    mirror.data[,c("x","y")] <- data.frame(x=tx,y=ty)
    I.dat <- .get_data(ribs[[1]], mirror.data)
    II.dat <- .get_data(ribs[[2]], mirror.data)

    pseudo.kokybe[i] <- dif.fun(I.dat, II.dat)

    if (n.m.keep){
      sim <- do.call(c,list(sim1,list(mirror.data)))
    }
  }
  if (n.m.keep){
    return(list(pseudo.kokybe,sim))
  } else {
    return(list(pseudo.kokybe))
  }
}

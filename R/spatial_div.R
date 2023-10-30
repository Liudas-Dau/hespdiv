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
  perim_pts <- .perimeter_pts(polygon = rims[[testid]],n.pts = n.split.pts,
                              dst.pts = dst.pts)

  .visualise_splits.start(what = trace.object,
                          pnts.col, xy.dat, rims,perim_pts,
                          testid, study.pol.viz)

  #testavimui pjuviai paruosiami

  environment(.dif_fun) <- e
  if (S.rel.crit > 0){
    S_org <- abs(pracma::polyarea(x=rims[[testid]][,1],y=rims[[testid]][,2]))
  } else {
    S_org <- NULL
  }

  pairs_pts <- .pair_pts(perim_pts[[1]],polygon = perim_pts[[2]])
  maxdif <- c.Q.crit # first split minimum quality. P.crit
  any.split <- numeric()
  maxid <- 0
  cond <- TRUE
  if (nrow(pairs_pts) != 0){
    #pjaustymo ir testavimo ciklas
    for (i in 1:nrow(pairs_pts)){
      if (!cond) .visualise_splits.bad_straight(what = trace.object,
                                                level = trace.level,
                                                pairs_pts, message, i-1)
      cond <- TRUE
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
      id1 <- .get_ids(po, samp.xy)
      id2 <- .get_ids(virs, samp.xy)

      .visualise_splits.try_straight(what = trace.object, level = trace.level,
                                     pairs_pts, i)

      if (N.crit > 0){
        cond <- length(id1) > N.crit & length(id2) > N.crit
        if (!cond){
          message <- paste0("Not enough observations.",
                            "\nObtained: ", length(id1), ' and ', length(id2),
                            "\nRequired: >", N.crit)
          next
        }
      }

      if (N.rel.crit > 0){
        cond <- length(id1) / nrow(samp.xy) > N.rel.crit &
          length(id2) / nrow(samp.xy) > N.rel.crit
        if (!cond){
          message <- paste0("Too low proportion of observations.\nObtained: ",
                            round(length(id1) / nrow(samp.xy), 2),
                            ' and ', round(length(id2) / nrow(samp.xy), 2),
                            "\nRequired: >", N.rel.crit)
          next
        }
      }

      if (N.loc.crit > 0 | N.loc.rel.crit > 0){
        n.uni.1 <- nrow(unique(samp.xy[id1,]))
        n.uni.2 <- nrow(unique(samp.xy[id2,]))
        if (N.loc.crit > 0){
          cond <- n.uni.1 > N.loc.crit &
            n.uni.2 > N.loc.crit
          if (!cond){
            message <- paste0("Not enough locations.",
                              "\nObtained: ", n.uni.1,
                              ' and ', n.uni.2,
                              "\nRequired: >", N.loc.crit)
            next
          }
        }
        if (N.loc.rel.crit > 0){
          n.uni <- nrow(unique(samp.xy))
          cond <- n.uni.1 / n.uni > N.loc.rel.crit &
            n.uni.2 / n.uni > N.loc.rel.crit
          if (!cond){
            message <- paste0("Too low proportion of locations.",
                              "\nObtained: ", round(n.uni.1 / n.uni,2),
                              ' and ', round(n.uni.2 / n.uni, 2),
                              "\nRequired: >", N.loc.rel.crit)
            next
          }
        }
      }

      if (S.crit > 0 | S.rel.crit > 0 ){
        SpjuvioI <- abs(pracma::polyarea(x=virs[,1],y=virs[,2]))
        SpjuvioII <- abs(pracma::polyarea(x=po[,1],y=po[,2]))
        if (S.crit > 0) {
          cond <- SpjuvioI > S.cond & SpjuvioII > S.cond
          if (!cond){
            message <- paste0("One of the areas was too small.\n","Obtained: ",
                              round(SpjuvioI,2), ' and ', round(SpjuvioII,2),
                              "\nRequired: >", S.cond)
            next
          }
        }
        if (S.rel.crit > 0) {
          cond <- SpjuvioI / S_org  > S.rel.crit &
            SpjuvioII / S_org > S.rel.crit
          if (!cond){
            message <- paste0("Too low proportion of area.\n","Obtained: ",
                              round(SpjuvioI / S_org,2), ' and ',
                              round(SpjuvioII / S_org, 2),
                              "\nRequired: >", S.rel.crit)
            next
          }
        }
      }
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
        cond <- FALSE
      }
    }
    if (!cond) .visualise_splits.bad_straight(what = trace.object,
                                              level = trace.level,
                                              pairs_pts, message, i)
    assign(x = "str.split.quals" ,
           value = do.call(c,list(str.split.quals,list(any.split))),
           envir = e)

    if(length(any.split) > 0){

      mean.dif <- mean(any.split) # mean spatial heterogeneity
      sd.dif <- sd(any.split) # NA if 1 split
      # sd(any.split) -> anisotropy of heterogeneity
      # if all are equal, then complete randomness: no anisotropy - no splits present
      if( all(any.split==any.split[1]) & length(any.split) > 1  ){
        message(paste0(c(
          "All tested splits were of equal quality in polygon "
          , testid, ". A random split that meets subdivision criteria was",
          " selected as best.")))
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
      # skirtas sumazinti poligona pateikiama?
      del.id <- unique(unlist(pairs_pts[,6:7])) # visi prideti taskai (atrinkti, taip kad nesidubliuotu)
      dub.id <- which(duplicated(perim_pts[[2]][-nrow(perim_pts[[2]]),]))
      except <- unique(c(
        unlist(apply(pairs_pts[maxid,6:7],2,function(o) which(del.id %in% o))),
        which(del.id == 1), # kurie prideti taskai yra pirmas taskas arba geriausias
        which(apply(perim_pts[[2]][del.id,],1,function(o)
          any(o[1] == perim_pts[[2]][dub.id,1] &
                o[2] == perim_pts[[2]][dub.id,2]) ))
      )) # kurie taskai yra kartu ir dubliai su
      # poligono taskais.

      del.id <- del.id[-except] # pasalinam situos taskus is del.id

      if (length(dub.id) > 0){
        if(any(c(pairs_pts[maxid,6:7]) %in% dub.id)) {
          except.id <- which(dub.id %in% pairs_pts[maxid,6:7])
          dubs <- dub.id[-except.id]
          except.id <- dub.id[except.id]
          alter.ids <- which(perim_pts[[2]]$x.poly %in%
                               perim_pts[[2]]$x.poly[except.id] &
                               perim_pts[[2]]$y.poly %in%
                               perim_pts[[2]]$y.poly[except.id])
          alter.ids <- alter.ids[ !alter.ids %in%
                                    c(pairs_pts[maxid,6],pairs_pts[maxid,7])]
          dub.id <- c(dubs, alter.ids)
        }
      }
      del.id <- unique(c(del.id, dub.id))
      if (length(del.id) > 0){
        perim_pts[[2]] <- perim_pts[[2]][-del.id,] # istrinam del.id is poligono
        min.x.id <- pairs_pts[maxid,6] - sum(del.id<pairs_pts[maxid,6])
        max.x.id <- pairs_pts[maxid,7] - sum(del.id<pairs_pts[maxid,7]) # pakoreguojam min.x.id ir max.x.id pagal istrintu
        # tasku kieki
      } else {
        min.x.id <- pairs_pts[maxid,6]
        max.x.id <- pairs_pts[maxid,7]
      }

      if ( c.splits) {
        .visualise_splits.best_straight(what = trace.object,
                                        pairs_pts, maxid, maxdif)

        environment(.curvial_split) <- environment()
        best.curve <- .curvial_split(

          poly.x = perim_pts[[2]]$x.poly,
          poly.y = perim_pts[[2]]$y.poly,
          min.x.id = min.x.id,
          max.x.id = max.x.id,
          b = pairs_pts[maxid,5],
          samp.dat = samp.dat,
          c.X.knots = c.X.knots,
          c.Y.knots = c.Y.knots,
          N.crit = N.crit,
          N.rel.crit = N.rel.crit,
          N.loc.crit = N.loc.crit,
          N.loc.rel.crit = N.loc.rel.crit,
          S.cond = S.cond,
          S.rel.crit = S.rel.crit,
          S_org = S_org,
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
           poly.info[testid,"c.improv"] >= c.crit.improv){
          poly.info[testid,"is.curve"] <- TRUE
          maxdif <- best.curve[[2]]
        }
      }
      if (!.comp(maxdif, Q.crit) ){
        maxid <- 0
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
      up.pol <- .close_poly(
        open.poly =
          .split_poly(
            polygon = data.frame(x=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],
                                 y=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
            split_ids = c(min.x.id, max.x.id),
            min_id = 1,
            trivial_side = TRUE,
            poli_side = TRUE
          ),
        close.line = best.splitl
      )
      up.ids <- .get_ids(up.pol, samp.xy)
      up.xy <- .slicer.table(samp.xy, up.ids)
      up.dat <- .slicer(samp.dat, up.ids)

      do.pol <- .close_poly(
        open.poly =
          .split_poly(
            polygon = data.frame(x=perim_pts[[2]][,1][-nrow(perim_pts[[2]])],
                                 y=perim_pts[[2]][,2][-nrow(perim_pts[[2]])]),
            split_ids = c(min.x.id, max.x.id),
            min_id = 1,
            trivial_side = TRUE,
            poli_side = FALSE
          ),
        close.line = best.splitl
      )
      do.ids <- .get_ids(do.pol, samp.xy)
      do.xy <- .slicer.table(samp.xy, do.ids)
      do.dat <- .slicer(samp.dat, do.ids)

      assign(x = "rims" ,value = do.call(c,list(rims,list(up.pol))) ,envir = e)

      .spatial_div(up.dat,up.xy, root.id = testid)
      assign(x = "rims" ,value = do.call(c,list(rims,list(do.pol))) ,envir = e)
      .spatial_div(do.dat,do.xy, root.id = testid)


    } # else no adequate split
  } else {
    assign(x = "str.split.quals" ,
           value = do.call(c,list(str.split.quals,list(any.split))),
           envir = e)
    assign(x = "poly.info" ,value = rbind(poly.info, data.frame(
      root.id = root.id,
      n.splits = 0,
      n.obs = nrow(samp.xy),
      mean = NA, # NA if 0
      sd = NA, # NA if 1 or 0 split
      str.best = NA,
      str.z.score = NA, # NA if 1 or 0 split
      has.split = FALSE,
      is.curve = FALSE,
      crv.best = NA,
      crv.z.score = NA,
      c.improv = NA))
      , envir = e)
    if (!is.null(trace.level))
      cat(paste("No adequate pair of points on perimeter of the polygon.",
                "\nMaybe too irregular polygon, too detailed subdivisions",
                "or too low n.pts."))
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





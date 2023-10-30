# CODE USED IN THE STUDY TO IMPORT & FILTER DATA & PERFORM ANALYSES
library(readxl)
library(xlsx)
library(janitor)
library(spData)
library(sp)
library(hespdiv)
library(tidyverse)
# Importing raw paleodb data
# ALL data downloaded 2023 April 12
# You can find META data in data-raw directory (names end with "_meta.txt")
early_mio_r <- read.csv(".//data-raw//early_mio.txt", sep = ',')
mid_mio_r <- read.csv(".//data-raw//mid_mio.txt", sep = ',')
late_mio_r <- read.csv(".//data-raw//late_mio.txt", sep = ',')
EARLY_r <- read.csv(".//data-raw//EARLY.txt", sep = ',')
LATE_r <- read.csv(".//data-raw//LATE.txt", sep = ',')
MIOCENE_r <- read.csv(".//data-raw//MIOCENE.txt", sep = ',')

raw_l <- list(early_mio = early_mio_r,
              mid_mio = mid_mio_r,
              late_mio =late_mio_r,
              EARLY = EARLY_r,
              LATE = LATE_r,
              MIOCENE = MIOCENE_r)

# Are there any observations, that are not identified to at least species
# level?
lapply(raw_l, function(o){
  nrow(o[which(o$accepted_rank != 'species' &
                 o$accepted_rank != 'subspecies' ),])
})
# NO

# Do all occurrences have coordinates?
lapply(raw_l, function(o){
  any(is.na(o$lng)) | any(is.na(o$lat))
})
# YES

# checking precision of localization.
layout(matrix(1:6,nrow = 3))
for (i in  1:6){
  barplot(table(raw_l[[i]]$latlng_precision), las = 2,main= names(raw_l)[i])
}
# most are of ok precision (error up to 1 degress or lower), should filter some obs.

for (i in 1:6){
  # checking what were the basis of localization
  barplot(sort(prop.table(table(raw_l[[i]]$latlng_basis))), las =2, cex.names = 0.8,
          main= names(raw_l)[i])
  # most of the obs. were localized based on political unit. Since it is the
  # standard localization technique, it is OK.

}
prec_l <- lapply(raw_l,function(o) o[-which((o$latlng_precision %in%
                                               as.character(2:7))),])

# How much obs. were removed?
x <- sapply(1:6,function(i) {
  rbind(nrow(raw_l[[i]]),nrow(prec_l[[i]]),
        nrow(raw_l[[i]]) - nrow(prec_l[[i]]))
})

colnames(x) <- names(raw_l)
rownames(x) <- c("Raw","Prec", "Dif")
x # not too much

# Checking if filtering was succesful
for (i in  1:6){
  barplot(table(prec_l[[i]]$latlng_precision), las = 2,main= names(prec_l)[i])
} # Yes
for (i in 1:6){
  barplot(sort(prop.table(table(prec_l[[i]]$latlng_basis))), las =2, cex.names = 0.8,
          main= names(prec_l)[i])
  # most of the obs. are still localized based on political unit.
}


# filtering away marine or/and coastal animals.
filter.mams <- function(mammalia){
   mammalia %>%
    as_tibble() %>%
    clean_names() %>%
    # should filter marine mammals
    filter(!(order %in% c("Cetacea", "Sirenia", "Desmostyloidea","Desmostylia",
                          "Chiroptera")),
           !(family %in% c("Odobenida","Phocidae",
                           "Otariidae","Desmatophocidae", "Odobenidae",
                           "Allodesminae")),# marine/coastal Carnivores
           !(genus %in% c("Enaliarctos", "Pteronarctos", "Pacificotaria",
                          "Pinnarctidion")), # marine/coastal genus, but NO_FAMILY_SPECIFIED

    )
}
ter_l <- lapply(prec_l, filter.mams)

# How much obs. left?
y <- sapply(1:6,function(i) {
  rbind(nrow(ter_l[[i]]),
        nrow(prec_l[[i]]) - nrow(ter_l[[i]]),
        nrow(raw_l[[i]]) - nrow(ter_l[[i]]))
})
colnames(y) <- names(raw_l)
rownames(y) <- c("Ter","Ter_Pre", "Ter_Raw")
x <- rbind(x,y)
x
# Enough.
# Are there any taxon names, containing at least 2 white spaces (or 3 words)?
suspicious_names <- lapply(ter_l, function(o) {
  unique((o %>% filter(str_detect(identified_name,
                                  pattern = '[:space:].*[:space:]')))$identified_name)})

suspicious_names[[1]]
suspicious_names[[2]]
suspicious_names[[3]]
suspicious_names[[4]]
suspicious_names[[5]]
suspicious_names[[6]]
# Yes there are, but it seems that most of them are new species or synonyms

# Are there any mentions of taxon names containing at least 2 white spaces (or
# 3 words) that do not have 'n. sp.' in it ("n. gen." not checked, because if
# taxon belongs to new genus it also should be a 'm. sp.')?
not_new_species <- lapply(ter_l, function(o) unique((o %>% filter(str_detect(identified_name,
                                                                             pattern = '[:space:].*[:space:]') &
                                                                    !str_detect(identified_name,
                                                                                pattern = 'n\\. sp\\.')))$identified_name))
not_new_species
# Yes there are, and all of them are synonym cases


# checking orders
lapply(ter_l, function(o) unique(o$order))

# are all species left are terrestrial?
lapply(ter_l, function(o) unique(o$taxon_environment))
# YES.

# Importing contiguous US coordinate data that will be used as the study
# polygon
# DO NOT RUN (as data depends on "spData" package we save it independently):
#data("world") # run 2023 April 13
# US polygon coordinates data.frame
#US <- world$geom[[5]][[10]][[1]]
# visual validation
#layout(1)
#plot(world$geom[[5]][[10]][[1]],type = 'l',xlab = 'x', ylab = "y")
# creating data.frame with coordinates of study area polygon:
#us_imp <- data.frame(x = US[,1], y =US[,2])
# Saving US polygon:
#write.xlsx(x = us_imp,file =  ".//data-raw//us.xlsx", col.names = TRUE,
#sheetName = "US polygon", append = FALSE,row.names = FALSE)
# Importing:
us <- data.frame(
  read_excel(".//data-raw//us.xlsx", sheet = "US polygon", col_names = TRUE)
)
# was the save successful?
#identical(us,us_imp)
# yes

# FIltering obs. by us polygon
us_l <- lapply(ter_l, function(o) {
  o[which(sp::point.in.polygon(pol.x = us[,1],
                                                  pol.y = us[,2],
                                                  point.x = o$lng,
                                                  point.y = o$lat) != 0),]
})
# Checking if there were any outliers
z <- sapply(1:6,function(i) {
  rbind(nrow(us_l[[i]]),
    nrow(ter_l[[i]]) - nrow(us_l[[i]]),
        nrow(prec_l[[i]]) - nrow(us_l[[i]]),
        nrow(raw_l[[i]]) - nrow(us_l[[i]]))
})
colnames(z) <- names(raw_l)
rownames(z) <- c("Us","Ter_US", "Pre_US","US_Raw")
x <- rbind(x,z)
x
# Yes, there were 8 obs from early miocene
# Where are they?
layout(1)
plot(us,type = 'l')
points(ter_l[[1]][which(sp::point.in.polygon(pol.x = us[,1],
                                             pol.y = us[,2],
                                             point.x = ter_l[[1]]$lng,
                                             point.y = ter_l[[1]]$lat) == 0),
                  c("lng","lat")],
       pch=  19, col=2 )
# Florida coast, as there only 8 obs., we ignore them.

# Are there any geogcomments that mentions uncertain?
lapply(us_l,function(o) (o %>% filter(str_detect(geogcomments,
                                     pattern = '.*uncertain.*')))$geogcomments)
# some obs. mentions uncertain county assignment. As counties are relatively
# small political units considering the scale of the analysis and purpose of
# this analysis, these observations are not removed.

# Because "major" obs. attribution to time bins rule was used when importing data
# from Paleobiology database, observations that were broadly dated "Miocene",
# will be assigned to a longer time interval (EARLY MIOCENE) in case when
# Miocene is divided into two parts. As these obs. do not have precise dating,
# they should be removed from analysis.

lapply(us_l, function(o) any("Miocene" %in% o$early_interval))
# These broadly dated obs. are ok for all Miocene (6th element of us_l)
# But not ok for Early Miocene.
# Here you can see that there are > 500 obs. dated "Miocene" in Early Miocene
# dataset, due to the use of major rule.
barplot(table(us_l[[4]]$early_interval),las =2,cex.names = .75)

# E.g. Mammuts evolved in Late Miocene, but due to the "major" rule they were
# assigned to Early Miocene:
us_l[[4]]$genus[1341]

#  Just to make sure that there any other surprises left:
identical(which((us_l[[4]]$max_ma -  us_l[[4]]$min_ma) ==
                   max(us_l[[4]]$max_ma -  us_l[[4]]$min_ma)),
           which(us_l[[4]]$early_interval == "Miocene"))

# So let's remove all obs. in Early Miocene dated "Miocene"

us_l[[4]] <- us_l[[4]][-which(us_l[[4]]$early_interval == "Miocene"),]
x[7,4] - nrow(us_l[[4]])
# 574 obs. removed.

# Let's see what families dominate in the filtered data sets.
layout(1:3)
for(i in 1:6 ){
  if (i <= 3){
    limy <- 500
  } else {
    if (i <= 5){
      if (i ==4) layout(1:2)
      limy <- 600
    } else {
      layout(1)
      limy <- 1100
    }

  }
  barplot(sort(table(us_l[[i]]$family)),ylim = c(0,limy),las = 2,
          cex.names = 0.75, main = names(us_l)[i])
  text(x = seq(0.5 ,length(unique(us_l[[i]]$family))*1.19,
               length.out = length(unique(us_l[[i]]$family))),
       y = sort(table(us_l[[i]]$family))+10,
       labels = sort(table(us_l[[i]]$family)),
       srt = 90,adj =c(0,0.5),cex = 0.7)
}

#Visually checking data
layout(cbind(1:3,4:6))
for ( i in 1:6){
  plot(us,type='l', main = colnames(x)[i])
  points(us_l[[i]]$lng,us_l[[i]]$lat,col=2)
}
# the main regions have obs
# missing late Miocene obs in North East coast

# Downloading and saving references of the occurrences that will be used in the
# study.
# DO NOT RUN:
# url2 <- "https://paleobiodb.org/data1.2/occs/refs.csv?base_name=Mammalia&rank=species&pres=regular&interval=Miocene&cc=US&state=!Alaska&select=occs"
# refs <- as.data.frame(read_csv(file = url2)) # downloaded 2023 April 13
# # Finding references of filtered occurrences in refs by matching
# # reference numbers:
# match_id <- lapply(us_l,function(o) match(o$reference_no,
#                       refs$reference_no))
# validating:
# sapply(1:6, function(i) all.equal(refs$reference_no[match_id[[i]]],
#           us_l[[i]]$reference_no)) # TRUE
#
# # filtering all refs by match_id:
# f_revs <- lapply(match_id,function(o) refs[o,])
# # saving:
# for (i in 1:6){
# write.xlsx(x = f_revs[[i]], file = paste0(".//data-raw//occ_references_",i,".xlsx"),
#                                           col.names = TRUE,
#                                           sheetName = "Matched with occurences",
#                                           append = FALSE,row.names = FALSE)
# # saving unique references:
# write.xlsx(x = unique(f_revs[[i]]), file = paste0(".//data-raw//occ_references_",i,".xlsx"),
#            col.names = TRUE,row.names = FALSE,
#            append = TRUE, sheetName = "Unique")
# }

# Saving workspace:
# DO NOT RUN:
# save.image("####") study_data.RData
#
# mio_mams <- us_l[[6]]
# There are some NonASCII
# apply(mio_mams,2,function(x){sum(grepl("[^ -~]",x))})[
#  which(apply(mio_mams,2,function(x){any(grepl("[^ -~]",x))}))]
# mio2 <- mio_mams
# for ( i in 1: ncol(mio_mams))
#  for (r in 1:nrow(mio_mams))
#    if(grepl("[^ -~]",mio_mams[r,i]))
#      mio2[r,i] <- gsub('[^\x20-\x7E]', '', mio_mams[r,i])
# all.equal(mio_mams,mio2) # NonASCII characters successfully removed
# mio_mams <- mio2
# usethis::use_data(mio_mams)
# usethis::use_data(us)


# species <- mio_mams$accepted_name # taxa names
# sp_coords <- data.frame(x = mio_mams$lng, y = mio_mams$lat)
# hd <- hespdiv(data = species, xy.dat = sp_coords, study.pol = us)
# usethis::use_data(hd)


# set.seed(2) # seed is used to obtained the same result of an experiment with random properties.
# hsa_res <- hsa(obj = hd,
#               n.runs = 100, # 100 alternative hespdiv re-reruns
#               n.split.pts = 8:30, # number of split-points determines fit to data of straight split-lines
#               same.n.split = FALSE, # split-point placement regularity determines whether scale of analysis changes depending on order of subdivision
#               c.splits = FALSE,  # Argument controls whether curves are generated in attempt to increase the performance of the best linear split-lines
#               c.X.knots = 3:8, # Controls the number of wiggles in generated curves, determining their fit to data
#               c.Y.knots = 5:15, # Controls the number of different shapes each curve wiggle can achieve, also determining their fit to data
#               c.fast.optim = TRUE, # Determines the optimization algorithm of non-linear split-lines
#               use.chull = FALSE) # Determines whether the convex polygon of occurences is used as a study area polygon. If not, study area polygon becomes the provided US polygon.
# usethis::use_data(hsa_res)

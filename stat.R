library(dplyr)
galkposrel=769500
galkposhs=816000
rel_length=4629891
hs_length=4634538

#diutsibutinb recombinaison

nombre_recombinaison <- function(recombinant){
  return(table(recombinant))
}


nombre_recombinaison(Fragment_position$recombinant) > 1

#distance_length
Fragment_position=read.csv("fragment_pos")
merge_frag <- function(Fragment_position, recomb){
  tmp <- Fragment_position[Fragment_position$recombinant==recomb,]
  Fragment_position <- Fragment_position[Fragment_position$recombinant!=recomb,]
  tmp[tmp$end==min(tmp$end),]$length <- tmp[tmp$begin==max(tmp$begin),]$length + tmp[tmp$end==min(tmp$end),]$length
  max_begin <- tmp[tmp$begin==max(tmp$begin),]$begin
  tmp <- tmp[tmp$begin!=max(tmp$begin),]
  tmp[tmp$end==min(tmp$end),]$begin <- max_begin
  return(rbind(Fragment_position, tmp))
}

Fragment_position<-merge_frag(Fragment_position, 18)
Fragment_position<-merge_frag(Fragment_position, 47)


distance_from_gal <- function(Fragment_position, galkpos, length_recipient){
  Fragment_position[Fragment_position$begin > (rel_length/2),]$begin = Fragment_position[Fragment_position$begin > (rel_length/2),]$begin-length_recipient
  Fragment_position[Fragment_position$end > (rel_length/2),]$end = Fragment_position[Fragment_position$end > (rel_length/2),]$end-length_recipient
  bloc_gal <- Fragment_position[(Fragment_position$begin < galkpos & Fragment_position$end > galkpos),]
  periph_gal <- Fragment_position[(Fragment_position$end < galkpos | Fragment_position$begin > galkpos),]
  
  periph_gal <- merge(periph_gal, select(bloc_gal, recombinant, begin), by="recombinant")
  periph_gal <- merge(periph_gal, select(bloc_gal, recombinant, end), by="recombinant")
  
  periph_gal$distance_1 = abs(periph_gal$begin.y - periph_gal$end.x)
  periph_gal$distance_2 = abs(periph_gal$end.y - periph_gal$begin.x)
  periph_gal$distance=apply(X=subset(periph_gal, select=c(distance_1, distance_2)), MARGIN=1, FUN=min)
  bloc_gal$distance = 0
  periph_gal=rename(periph_gal, begin=begin.x)
  periph_gal=rename(periph_gal, end=end.x)
  periph_gal=subset(periph_gal, select=-c(begin.y, end.y, distance_1, distance_2))
  
  return(rbind(bloc_gal, periph_gal))
}


Fragment_position = rbind(distance_from_gal(Fragment_position[Fragment_position$recombinant < 29,], galkposrel, rel_length),
      distance_from_gal(Fragment_position[Fragment_position$recombinant >= 29,], galkposhs, hs_length))

plot(Fragment_position$distance, Fragment_position$length, xlim=c(0, 3e6))
dist_nb_recom = nombre_recombinaison(Fragment_position$recombinant)
#couververture gal



couverture_gal <- function(Fragment_position, galkpos){
  Fragment_position <- Fragment_position[(Fragment_position$begin < galkpos & Fragment_position$end > galkpos),]
  range = max(Fragment_position$end) - min(Fragment_position$begin) + 1
  Donor<-rep(0, range)
  for(recombinant in unique(Fragment_position$recombinant)){
    begin = Fragment_position[Fragment_position$recombinant==recombinant,]$begin
    end = Fragment_position[Fragment_position$recombinant==recombinant,]$end
    Donor[(begin-min(Fragment_position$begin)+1):(end-min(Fragment_position$begin)+1)] = Donor[(begin-min(Fragment_position$begin)+1):(end-min(Fragment_position$begin)+1)] + 1
  }
  return(Donor)
}


donor_rel = couverture_gal(Fragment_position[Fragment_position$recombinant < 29,], galkposrel)
donor_hs = couverture_gal(Fragment_position[Fragment_position$recombinant >= 29,], galkposhs)

plot(1:length(donor_hs), donor_hs, ylim=c(0, 20))
plot(1:length(donor_rel), donor_rel, ylim=c(0, 20))
abline(v=galkposhs-min(bloc_gal$begin)+1)



hist(Fragment_position$length, breaks= 15, xlab = "Longueur des blocs",
     main="Distribution de la longueur des blocs recombinés")


# inter <- data.frame(rec = numeric(13), long = numeric(13))
# inter$rec[1] <- Fragment_position$recombinant[1]
# inter$long[1] <- NA
# for(i in 2:14){
#   if (Fragment_position$recombinant[i] == Fragment_position$recombinant[i-1]){
#     inter$rec[i-1] <- Fragment_position$recombinant[i]
#     inter$long[i-1] <- Fragment_position$begin[i]- Fragment_position$end[i-1]}
#   else{
#     inter$rec[i-1] <- Fragment_position$recombinant[i] -1
#     inter$long[i-1] <- NA}
# }
#
# print(na.omit(inter))
# hist(inter$long,breaks=100, xlab = 'bp', main = "Distance entre les blocs recombinés")


table_files = read.table("position_files.txt")
data.frame(K12=table_files$V1[1:(length(table_files$V1)/2)], REL606=table_files$V1[((length(table_files$V1)/2)+1):(length(table_files$V1))])






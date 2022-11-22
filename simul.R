galkposrel=769500
galkposhs=816000
rel_length=4629891
hs_length=4634538

simulation_distance <- function(Fragment_position, effectif, galkpos, length){
  bloc_gal <- Fragment_position[(Fragment_position$begin < galkpos & Fragment_position$end > galkpos),]
  periph_block <- Fragment_position[!(row.names(Fragment_position) %in% row.names(bloc_gal)),]
  distance=c()
  for(i in 1:effectif){
    recombinant = sample(unique(periph_block$recombinant), 1)
    begin = bloc_gal[bloc_gal$recombinant==recombinant,]$begin
    end = bloc_gal[bloc_gal$recombinant==recombinant,]$end
    vec = c(-(length/2):begin,end:(length/2))
    positions_tmp = sample(vec, nrow(periph_block[periph_block$recombinant == recombinant,]))
    print(positions_tmp)
    distance = c(distance, abs(positions_tmp[positions_tmp < begin] - begin), abs(positions_tmp[positions_tmp > end] - end))
  }
  return(distance)
}


simulation_longueur <- function(Fragment_position, effectif, galkpos){
  bloc_gal <- Fragment_position[(Fragment_position$begin < galkpos & Fragment_position$end > galkpos),]
  periph_block <- Fragment_position[!(row.names(Fragment_position) %in% row.names(bloc_gal)),]
  bloc_gal <- bloc_gal[(bloc_gal$recombinant%in% periph_block$recombinant),]
  counter=0
  for(i in 1:effectif){
    recombinant = sample(bloc_gal$recombinant, 1)
    len_block = sample(Fragment_position[Fragment_position$recombinant==recombinant, ]$length, 1)
    if(len_block == max(Fragment_position[Fragment_position$recombinant==recombinant, ]$length)){
      counter = counter + 1
    }
  }
  return(counter/effectif)
}

vec = c((-gal_length/2):bloc_gal[bloc_gal$recombinant==12,]$begin,bloc_gal[bloc_gal$recombinant==12,]$end:(gal_length/2))
sample(vec, nrow(periph_block[periph_block$recombinant == 12,]))

print(simulation_distance(Fragment_position, 10))
simulation_longueur(Fragment_position, 50)

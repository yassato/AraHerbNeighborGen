coord = function(chr, pos) {
  if(length(pos)!=length(chr)) stop("chr and pos length differ")
  chr <- as.factor(chr)
  coord <- 0
  M <- 0
  for (i in 1:nlevels(chr)) {
    w <- (chr == levels(chr)[i])
    pos.c <- pos[w]
    coord[w] <- M + pos.c
    mx <- max(pos.c)
    M <- M + mx
  }
  coord <- coord/M
  return(coord)
}

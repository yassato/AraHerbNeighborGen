coord = function(chr, pos) {
  if(length(pos)!=length(chr)) stop("chr and pos length differ")
  chr <- as.factor(chr)
  coord <- 0
  M <- 0
  tic <- numeric(nlevels(chr))
  for (i in 1:nlevels(chr)) {
            w <- (chr == levels(chr)[i])
            pos.c <- pos[w]
            coord[w] <- M + pos.c
            mx <- max(pos.c)
            tic[i] <- M + mx/2
            M <- M + mx
            }
  coord <- coord/M
  tic <- tic/M
  return(list(coord=coord,tic=tic))
}

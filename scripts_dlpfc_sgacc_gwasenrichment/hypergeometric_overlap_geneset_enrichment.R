# commands to test hypergeometric enrichment 
# between DGEs in psychiatric disorders and DGEs between DLPFC and sgACC 

getRF <- function(list1, list2 ) {
  # http://nemates.org/MA/progs/representation.stats.html 
  both <- length(unique(intersect(list1, list2)))
  g1 = length(unique(list1))
  g2 = length(unique(list2))
  
  pval <- phyper(both-1, g1, totalgenome-g1, g2, lower.tail = FALSE, log.p = FALSE)
  rf <- (both-1)/((g1*g2)/totalgenome)
  return(rf)
  
}
getP <- function(list1, list2 ) {
  # http://nemates.org/MA/progs/representation.stats.html 
  both <- length(unique(intersect(list1, list2)))
  g1 = length(unique(list1))
  g2 = length(unique(list2))
  
  pval <- phyper(both-1, g1, totalgenome-g1, g2, lower.tail = FALSE, log.p = FALSE)
  rf <- both/((g1*g2)/totalgenome)
  return(pval)
  
}
getInter <- function(list1, list2 ) {
  # http://nemates.org/MA/progs/representation.stats.html 
  both <- length(unique(intersect(list1, list2)))
  g1 = length(unique(list1))
  g2 = length(unique(list2))
  
  #pval <- phyper(both-1, g1, totalgenome-g1, g2, lower.tail = FALSE, log.p = FALSE)
  #rf <- both/((g1*g2)/totalgenome)
  return(both)
  
}
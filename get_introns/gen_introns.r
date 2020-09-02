suppressPackageStartupMessages(library(GenomicFeatures))

txdb <- makeTxDbFromGFF('hg19.ncbiRefSeq.gtf')

all.introns <- intronicParts(txdb)

df <- data.frame(seqnames=seqnames(all.introns),
  starts=start(all.introns)-1,
  ends=end(all.introns),
  names=c(rep(".", length(all.introns))),
  scores=c(rep("1", length(all.introns))),
  strands=strand(all.introns))

write.table(df, file="all_introns.bed", quote=F, sep="\t", row.names=F, col.names=F)

' run_sequenza.R

Usage: 
run_sequenza.R -i INPUT -p PREFIX -o OUTPUT

Options:
-i --INPUT input      input .seqz file
-p --PREFIX prefix    sample prefix
-o --OUTPUT output    output directory
' -> doc

library(docopt)
library(dplyr)
library(sequenza)
args <- docopt(doc)
infile=args$INPUT
outfile=args$OUTPUT
prefix=args$PREFIX
seqz.data <- sequenza.extract(infile, chromosome.list=c(1:22, "X"))
#seqz.data <- sequenza.extract(infile, chromosome.list="22")
#gc.stats <- gc.norm(x = seqz.data$depth.ratio,gc = seqz.data$GC.percent)
#gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
#seqz.data$adjusted.ratio <- seqz.data$depth.ratio/gc.vect[as.character(seqz.data$GC.percent)]

fit=sequenza.fit(seqz.data)
cint=get.ci(fit)

sequenza.results(sequenza.extract = seqz.data, cp.table = fit, sample.id = prefix, out.dir=outfile)

#result=data.frame("ploidy" = cint$max.ploidy, "tc" = cint$max.cellularity)

#write.table(result, outfile, sep = "\t", row.names = F)
context('event detection functions')
colo829 <- readVcf(system.file("extdata", "diploidSV.vcf", package = "RTDetect"))

#RT detection
#vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#gr <- breakpointRanges(vcf, nominalPosition=TRUE)
#rt <- rtDetect(gr, genes, maxgap=30, minscore=0.6)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes <- TxDb.Hsapiens.UCSC.hg19.knownGene

test_that("RT detection returns two Granges", {
    gr <- StructuralVariantAnnotation::breakpointRanges(colo829, nominalPosition=TRUE)
    rt <- rtDetect(gr, genes, maxgap=50, minscore=0.3)
    if (length(rt)>0) {
        expect_equal(rep(2, length(rt)), unname(sapply(rt, length)))
    }
    
})


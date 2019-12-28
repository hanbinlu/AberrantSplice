package main

import (
	//"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
	"fmt"
	"github.com/Hanbin/AberrantSplice/Internal/gtfparser"
	"os"
	"time"
	//"sync"
)

func main() {
	gtf := os.Args[1]
	//geneset := []string{"sTCF3"}
	//geneset := genodatastruct.Allgene
	//start := time.Now()
	//result := gtfparser.Parsegtf(gtf, geneset)
	//fmt.Printf("Sequential parse takes %v\n%v\n", time.Since(start), result["ENSG00000071564.16"])
	start := time.Now()
	genes := gtfparser.ParsegtfConcurrent(gtf)
	fmt.Printf("Concurrent parse %v genes takes %v\n", len(genes), time.Since(start))

	//introns := make(map[string][]genodatastruct.Coor)
	//for _, gene := range genes {
	//	introns[gene.Chromosome] = append(introns[gene.Chromosome], gene.IntervalOfExons()...)
	//}
	//merge overlaping introns and sort output
	//for chro := range introns {
	//	regions := introns[chro]
	//	sort.Slice(regions, func(i, j int) bool { return regions[i].Start <= regions[j].Start })
	//	if len(regions) > 0 {
	//		fmt.Println("%v: %v", chro, regions[0:9])
	//	} else {
	//		log.Println("%v has no introns", chro)
	//	}
	//}
}

//func AssignMappedReads(cigar *genodatastruct.CIGAR, genes map[string]*genodatastruct.Gene) {
//
//}

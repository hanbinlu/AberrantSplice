package main

import (
	"os"
	"sync"

	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
	"github.com/Hanbin/AberrantSplice/Internal/gtfparser"
	"github.com/Hanbin/AberrantSplice/Internal/samparser"
	"github.com/Hanbin/AberrantSplice/scripts/splicetype"
)

func main() {
	gtf, sam := os.Args[1], os.Args[2]
	genes := gtfparser.ParsegtfConcurrent(gtf)
	index := splicetype.SortGeneMap(genes)
	mapqfilter := func(s genodatastruct.SamRec) bool {
		if s.MAPQ > 30 {
			return true
		}
		return false
	}
	samchan := samparser.ParseSam(sam, mapqfilter)
	normal, intronInc, exonSkip := 0, 0, 0
	var out []chan string
	nworker := 6
	for i := 0; i < nworker; i++ {
		o := make(chan string)
		out = append(out, o)
		worker := splicetype.RMTConstructor{
			In:    samchan,
			Out:   o,
			Index: index,
			Genes: genes,
		}
		go worker.Construct()
	}
	//merge chan
	var wg sync.WaitGroup
	wg.Add(nworker)
	mergechan := make(chan string)
	output := func(c chan string) {
		for v := range c {
			mergechan <- v
		}
		wg.Done()
	}
	for _, o := range out {
		go output(o)
	}
	go func() {
		wg.Wait()
		close(mergechan)
	}()
	//take results
	for s := range mergechan {
		if s == "normal" {
			normal++
		}
		if s == "intronInclusion" {
			intronInc++
		}
		if s == "exonSkipping" {
			exonSkip++
		}
		//if s == "fatal" {
		//	log.Fatalln("Exception read")
		//}
	}

	println("Normal reads #", normal)
	println("Intron Inclusion reads #", intronInc)
	println("Exon skipping reads #", exonSkip)
}

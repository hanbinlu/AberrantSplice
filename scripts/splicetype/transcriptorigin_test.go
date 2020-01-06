package splicetype

import (
	"os"
	"testing"

	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
	"github.com/Hanbin/AberrantSplice/Internal/gtfparser"
	"github.com/Hanbin/AberrantSplice/Internal/samparser"
)

func TestTranscriptOrigin(t *testing.T) {
	gtf := os.Args[1]
	sam := os.Args[2]
	stabletest(gtf, sam)
	//stat := 0
	//lastresult := []string{}
	//for i := 0; i < 6; i++ {
	//	result := stabletest(gtf, sam)
	//	if stat != 0 {
	//		compareresult(result, lastresult)
	//	}
	//	if stat != 0 && len(result) != stat {
	//		log.Fatalln("Break")
	//	}
	//	stat = len(result)
	//	lastresult = result
	//}
}

func stabletest(gtf, sam string) {
	//result := []string{}
	genes := gtfparser.ParsegtfConcurrent(gtf)
	index := SortGeneMap(genes)
	mapqfilter := func(s genodatastruct.SamRec) bool {
		if s.MAPQ > 30 {
			return true
		}
		return false
	}
	samchan := samparser.ParseSam(sam, mapqfilter)
	total, cnt := 0, 0
	var out []chan []int
	for i := 0; i < 6; i++ {
		o := make(chan []int)
		out = append(out, o)
		worker := RMTConstructor{
			in:    samchan,
			out:   o,
			index: index,
			genes: genes,
		}
		go worker.Construct()
	}
	for _, o := range out {
		raw := <-o
		cnt += raw[0]
		total += raw[1]
	}
	//for samrec := range samchan {
	//	mr := ReadMapTranscriptome{
	//		Chromosome: samrec.Chromosome,
	//		Segment:    samrec.RegionAligned(),
	//	}
	//	mr.InvolvedGeneLoci(index)
	//	if len(mr.GeneLoci) > 1 {
	//		cnt++
	//	}
	//}
	println(cnt, "######")
}

//func compareresult(A, B []string) {
//	if len(A) == len(B) {
//		for i, v := range A {
//			if v != B[i] {
//				log.Fatalln("######Unmatch result", v, " with ", B[i])
//			}
//		}
//		println("Exactly the same")
//	} else if len(A) > len(B) {
//		j := 0
//		for i, v := range A {
//			if v != B[j] {
//				log.Print("Unmatch result ", v, " with ", B[i], "Resolved ")
//				if B[j] == A[i+1] {
//					log.Println("It is an insertion")
//					continue
//				}
//			}
//			j++
//		}
//	} else {
//		compareresult(B, A)
//	}
//}

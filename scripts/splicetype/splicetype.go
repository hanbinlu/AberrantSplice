package splicetype

import (
	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
)

//import (
//	"math"
//)
//"sort"

//"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"

type TranOri string

const Exonic, Intronic, Crossonic TranOri = "E", "I", "C"

//Segment is exonic, intronic or crossonic (spanning intron and exon)
func (m *TranCoor) ClassSeg() TranOri {
	if len(m.IntronID) == 0 && len(m.ExonID) != 0 {
		return Exonic
	} else if len(m.IntronID) != 0 && len(m.ExonID) == 0 {
		return Intronic
	} else {
		return Crossonic
	}
}

//Values: normal, intronic, intronInclusion, exonSkipping, truncExon...
type SpliceType string

func (mr *ReadMapTranscriptome) SpliceType() string {
	if len(mr.Segment) == 1 { //no junction
		sumflag := Crossonic //lowest value of TranOri type
		for _, trancoor := range mr.MapTran {
			flag := trancoor[0].ClassSeg()
			if flag == Exonic {
				sumflag = Exonic
				break
			} else if flag == Intronic {
				sumflag = Intronic
			} //flag == Crossonic no change the sumflag value
		}
		if sumflag == Exonic {
			return "normal"
		} else if sumflag == Intronic {
			return "intronic"
		} else {
			return "intronInclusion"
		}
	} else if len(mr.Segment) > 1 { // splited read
		classify := []string{}
		for _, trancoors := range mr.MapTran {
			temp := []TranOri{}
			for _, tc := range trancoors {
				temp = append(temp, tc.ClassSeg())
			}
			if len(temp) == 0 {
				continue
			}
			//detect intron inclusion and exon skipping
			if All(temp, func(tr TranOri) bool { return tr == Exonic }) {
				tag := "normal"
				for i := 0; i < len(trancoors)-1; i++ {
					if trancoors[i+1].ExonID[0]-trancoors[i].ExonID[0] > 1 {
						tag = "exonSkipping"
						break
					}
				}
				classify = append(classify, tag)
			} else {
				classify = append(classify, "intronInclusion")
			}
		}
		if AnyString(classify, func(s string) bool { return s == "normal" }) {
			return "normal"
		} else if AnyString(classify, func(s string) bool { return s == "exonSkipping" }) {
			return "exonSkipping"
		} else if AnyString(classify, func(s string) bool { return s == "intronInclusion" }) {
			return "intronInclusion"
		}
	}
	return "No Class"
}

//Goroutine infrastruture to generate
type RMTConstructor struct {
	In    <-chan genodatastruct.SamRec
	Out   chan string
	Genes map[string]*genodatastruct.Gene
	Index map[string]*GeneMapIndex
}

func (w *RMTConstructor) Construct() {
	//total := 0
	for samrec := range w.In {
		mr := ReadMapTranscriptome{
			Chromosome: samrec.Chromosome,
			Strand:     samrec.Strand(),
			Segment:    samrec.RegionAligned(),
		}
		mr.InvolvedGeneLoci(w.Index)
		mr.MapToTran(w.Genes)
		if len(mr.MapTran) == 0 {
			continue
		}
		//total++
		w.Out <- mr.SpliceType()
	}
	//println("I have processed ", total)
	close(w.Out)
}

func AnyString(vs []string, f func(string) bool) bool {
	for _, v := range vs {
		if f(v) {
			return true
		}
	}
	return false
}
func Any(vs []TranOri, f func(TranOri) bool) bool {
	for _, v := range vs {
		if f(v) {
			return true
		}
	}
	return false
}
func All(vs []TranOri, f func(TranOri) bool) bool {
	for _, v := range vs {
		if !f(v) {
			return false
		}
	}
	return true
}

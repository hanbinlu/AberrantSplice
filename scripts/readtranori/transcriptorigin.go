package readtranori

import (
	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
	"log"
	"sort"
)

//Wrapper for gene map index for each chromosome
//non-overlapping loci merge intersecting gene loci
//for fast locate searching region of reads
type GeneMapIndex struct {
	GeneLoci       []GeneCoorPair
	NonOverlapLoci []genodatastruct.Coor
}
type GeneCoorPair struct {
	GeneID string
	Locus  genodatastruct.Coor
}

//Organize gene map by chromasome and sort
func SortGeneMap(Genes map[string]*genodatastruct.Gene) map[string]*GeneMapIndex {
	GMI := map[string]*GeneMapIndex{}
	for geneid, gene := range Genes {
		if _, ok := GMI[gene.Chromosome]; !ok {
			GMI[gene.Chromosome] = &GeneMapIndex{}
		}
		GMI[gene.Chromosome].GeneLoci = append(GMI[gene.Chromosome].GeneLoci, GeneCoorPair{geneid, gene.Coordinate})
	}
	for _, gmi := range GMI {
		//sort the gene loci by start coor
		sort.SliceStable(gmi.GeneLoci, func(i, j int) bool {
			return gmi.GeneLoci[i].Locus.Start <= gmi.GeneLoci[j].Locus.Start
		})
		//compress into non overlapping loci
		tomerge := []genodatastruct.Coor{}
		for _, v := range gmi.GeneLoci {
			tomerge = append(tomerge, v.Locus)
		}
		gmi.NonOverlapLoci = genodatastruct.MergeRegions(tomerge)
	}
	return GMI
}

//Map mapped read to the transcriptomic origin (RMT)
type ReadMapTranscriptome struct {
	Chromosome string
	Segment    []genodatastruct.Coor
	GeneLoci   []string
	ExonOri    map[string][]int //transcriptname->matched exon number of each segment
}

//Searching for gene loci that Intersect with any of the segment
func (mr *ReadMapTranscriptome) InvolvedGeneLoci(index map[string]*GeneMapIndex) {
	geneset := map[string]bool{}
	//check index
	if _, ok := index[mr.Chromosome]; !ok {
		mr.GeneLoci = append(mr.GeneLoci, "")
		return
	}

	//Narrow down the searching range by binary searching sorted index
	end := sort.Search(len(index[mr.Chromosome].GeneLoci), func(i int) bool {
		return index[mr.Chromosome].GeneLoci[i].Locus.Start > mr.Segment[len(mr.Segment)-1].End
	})
	if end == len(index[mr.Chromosome].GeneLoci) {
		end--
	}
	//searching from the non-overlapping loci
	boundaryIdx := sort.Search(len(index[mr.Chromosome].NonOverlapLoci), func(i int) bool {
		return index[mr.Chromosome].NonOverlapLoci[i].End >= mr.Segment[0].Start
	}) - 1
	if boundaryIdx == -1 {
		boundaryIdx = 0
	}
	boundary := index[mr.Chromosome].NonOverlapLoci[boundaryIdx].End

	//searching backwards from the end gene loci
	for i := end; i > 0; i-- {
		genecoorpair := index[mr.Chromosome].GeneLoci[i]
		if genecoorpair.Locus.Start <= boundary {
			//the rest wont intersect with the read.
			//terminate searching
			break
		}
		for _, seg := range mr.Segment {
			if seg.Intersect(genecoorpair.Locus) {
				if _, ok := geneset[genecoorpair.GeneID]; !ok {
					//println(genecoorpair.GeneId)
					geneset[genecoorpair.GeneID] = true
					mr.GeneLoci = append(mr.GeneLoci, genecoorpair.GeneID)
				}
			}
		}
	}
	if len(mr.GeneLoci) == 0 { //no hit
		mr.GeneLoci = append(mr.GeneLoci, "")
	}
	return
}

//Matching the segment to exon
func (mr *ReadMapTranscriptome) LocateExon(genes map[string]*genodatastruct.Gene) {
	if mr.GeneLoci[0] == "" {
		return
	} else if len(mr.GeneLoci) == 0 {
		log.Fatalln("Call InvolvedGeneLoci() first")
	}
	exonori := map[string][]int{}
	for _, geneid := range mr.GeneLoci {
		//go over the transcripts of each related gene
		if _, ok := genes[geneid]; !ok {
			log.Fatalln(geneid, " is missing in GTF file")
		}
		transcripts := genes[geneid].Transcripts
		for _, t := range transcripts {
			//set the temp to be -1. -1 stands for no matching exon from this transcript
			temp := make([]int, len(mr.Segment))
			for x := range temp {
				temp[x] = -1
			}
			//compare each exon with each segement
			for i, exon := range t.Exons {
				for j, seg := range mr.Segment {
					if seg.Inside(exon) {
						temp[j] = i //i exon of the transcript contains segment j
						exonori[t.TranscriptName] = temp
					}
				}
			}
		}
	}
	mr.ExonOri = exonori
}

//Goroutine infrastruture to generate
type RMTConstructor struct {
	in    <-chan genodatastruct.SamRecPartial
	out   chan int
	genes map[string]*genodatastruct.Gene
	index map[string]*GeneMapIndex
}

func (w *RMTConstructor) Construct() {
	cnt := 0
	for samrec := range w.in {
		mr := ReadMapTranscriptome{
			Chromosome: samrec.Chromosome,
			Segment:    samrec.RegionAligned(),
		}
		mr.InvolvedGeneLoci(w.index)
		mr.LocateExon(w.genes)
		sum := 0
		for _, exonlist := range mr.ExonOri {
			for _, v := range exonlist {
				if v != -1 {
					sum++
				}
			}
		}
		if sum == 0 {
			cnt++
		}
	}
	w.out <- cnt
	close(w.out)
}

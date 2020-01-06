package splicetype

import (
	"log"
	"sort"

	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
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
	Strand     string
	Segment    []genodatastruct.Coor
	GeneLoci   []string
	MapTran    [][]TranCoor //transcriptname->matched exon number of each segment
}

//Searching for gene loci that Intersect with any of the segment
func (mr *ReadMapTranscriptome) InvolvedGeneLoci(index map[string]*GeneMapIndex) {
	geneset := map[string]bool{}
	//check index
	if _, ok := index[mr.Chromosome]; !ok {
		mr.GeneLoci = append(mr.GeneLoci, "No Chromosome")
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
			if seg.Inside(genecoorpair.Locus) { //only consider genic reads
				if _, ok := geneset[genecoorpair.GeneID]; !ok {
					//println(genecoorpair.GeneId)
					geneset[genecoorpair.GeneID] = true
					mr.GeneLoci = append(mr.GeneLoci, genecoorpair.GeneID)
				}
			}
		}
	}
	if len(mr.GeneLoci) == 0 { //no hit
		mr.GeneLoci = append(mr.GeneLoci, "Intergenic")
	}
	return
}

//Detail annotate the segment's Exon/Intron location in the transcriptome
type TranCoor struct {
	TranscriptName string
	IsIn           bool
	ExonID         []int
	IntronID       []int
}

//Determine a single segment's splice status
func (mr *ReadMapTranscriptome) MapToTran(genes map[string]*genodatastruct.Gene) {
	if mr.GeneLoci[0] == "No Chromosome" || mr.GeneLoci[0] == "Intergenic" {
		return
	} else if len(mr.GeneLoci) == 0 {
		log.Fatalln("Call InvolvedGeneLoci() first")
	}
	//tranloc := [][]TranCoor{}
	for _, geneid := range mr.GeneLoci {
		//go over the transcripts of each related gene
		if _, ok := genes[geneid]; !ok {
			log.Fatalln(geneid, " is missing in GTF file")
		}
		if mr.Strand != genes[geneid].Strand {
			continue
		}
		transcripts := genes[geneid].Transcripts
		for _, t := range transcripts {
			temp := make([]TranCoor, len(mr.Segment))
			//Initiate Temp
			emptyCnt := 0
			for i, seg := range mr.Segment {
				temp[i] = TranCoor{TranscriptName: t.TranscriptName}
				//which exon and intron
				temp[i].ExonID = t.WhichExonIntersect(seg)
				temp[i].IntronID = t.WhichIntronIntersect(seg)
				if len(temp[i].ExonID)+len(temp[i].IntronID) == 0 {
					emptyCnt++
				}
			}
			if emptyCnt == 0 {
				mr.MapTran = append(mr.MapTran, temp)
			}
		}
	}
}

//func (mr *ReadMapTranscriptome)

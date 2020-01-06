package gtfparser

import (
	"bufio"
	"log"
	"os"
	"strconv"
	"strings"
	"sync"

	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
)

//This scripts try to fatorize the work flow of parsing a gtf file to enable concurrent
//processing. The parsing can be splited into following computational units
//read files by line (I/O)
//split into fields and pre-split attributes
//process attributes into corresponding data structure
//update to the output data structure (Gene struct)
//To streamline the above steps we will create following goroutines
//read files by line, split into fields and gather all the records of a gene -->
//Process into a Gene struct --> Put into a map[string]*Gene

//goroutine to produce gene lines to process
type Field []string
type GatherGeneRecs struct {
	gtf string
	out chan []Field //splited fields of all the gtf line related to a gene
}

//TakeGeneLines streamming out gene records
func (w *GatherGeneRecs) TakeGeneLines() {
	go func() {
		gtfF, err := os.Open(w.gtf)
		if err != nil {
			log.Fatal(err)
		}
		defer gtfF.Close()

		scanner := bufio.NewScanner(gtfF)
		genelines := []Field{}
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, "#") { //skip comment lines
				continue
			}
			fields := Field(strings.Split(line, "\t"))
			//is gene?
			if FeatureType(fields[2]) == Geneline {
				if len(genelines) > 0 {
					w.out <- genelines
				}
				genelines = []Field{fields}
			} else {
				genelines = append(genelines, fields)
			}
		}
		close(w.out)
	}()
}

//goroutines to produce a Gene struct
type GeneMapUnit struct {
	gene   *genodatastruct.Gene
	geneid string
}
type MakeGeneStruct struct {
	recieve chan []Field
	out     chan GeneMapUnit
}

//
func (w *MakeGeneStruct) GenereateGene() {
	go func() {
		for lines := range w.recieve { //grab lines from a field-splited gene
			var gene *genodatastruct.Gene
			var geneid string
			for _, fields := range lines {
				attributes := parseAttr(fields[len(fields)-1])
				//transform gtf record into corresponding data structure, Gene, Transcript
				if FeatureType(fields[2]) == Geneline {
					gene = makeGene(fields, attributes)
					geneid = attributes["gene_id"]
				} else if FeatureType(fields[2]) == Transcriptline {
					gene.Transcripts = append(gene.Transcripts, makeTranscript(fields, attributes))
				} else if FeatureType(fields[2]) == Exonline {
					idx := len(gene.Transcripts) - 1
					start, _ := strconv.Atoi(fields[3])
					end, _ := strconv.Atoi(fields[4])
					gene.Transcripts[idx].Exons = append(gene.Transcripts[idx].Exons, genodatastruct.Coor{start, end})
				}
			}
			w.out <- GeneMapUnit{gene, geneid}
		}
		close(w.out)
	}()
}

func ParsegtfConcurrent(gtf string) map[string]*genodatastruct.Gene {
	//Digest into gene lines from gtf file for further process into Gene struct
	readgene := GatherGeneRecs{
		gtf: gtf,
		out: make(chan []Field, 100),
	}
	readgene.TakeGeneLines()
	//Recieve the digest lines and make gene struct
	//Spawn multiple workers to parse
	const n = 6
	genemakers := [n]MakeGeneStruct{}
	for i := 0; i < n; i++ {
		genemakers[i] = MakeGeneStruct{
			recieve: readgene.out,
			out:     make(chan GeneMapUnit),
		}
		genemakers[i].GenereateGene()
	}
	//Gather genes produced by genemakers and make final gene map
	Genes := make(map[string]*genodatastruct.Gene)
	mergechan := func() chan GeneMapUnit {
		var wg sync.WaitGroup
		out := make(chan GeneMapUnit)
		output := func(c chan GeneMapUnit) {
			for m := range c {
				out <- m
			}
			wg.Done()
		}
		wg.Add(n)
		for _, worker := range genemakers {
			go output(worker.out)
		}
		go func() {
			wg.Wait()
			close(out)
		}()
		return out
	}()
	for unit := range mergechan {
		//sort the exons by strand
		strand := unit.gene.Strand
		for _, transcript := range unit.gene.Transcripts {
			exons := transcript.Exons
			if strand == "+" {
				genodatastruct.SortCoors(exons, true)
			} else {
				genodatastruct.SortCoors(exons, false)
			}
			transcript.Introns = transcript.GenerateIntrons()
		}
		Genes[unit.geneid] = unit.gene
	}
	return Genes
}

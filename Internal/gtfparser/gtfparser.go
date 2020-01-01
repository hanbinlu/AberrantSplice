package gtfparser

import (
	"bufio"
	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
	"log"
	"os"
	"strconv"
	"strings"
)

//FeatureType is string variant denotate a line's feature content of gtf file
type FeatureType string

//Geneline, Trancriptline, Exonline are 3 major features recorded in Gtf
const Geneline, Transcriptline, Exonline FeatureType = "gene", "transcript", "exon"

//Parsegtf parses gtf file and return parsed gene in the geneset
//if geneset=["all"] (a constant) it will stores all the genes
func Parsegtf(gtf string, geneset []string) map[string]*genodatastruct.Gene {
	gtfF, err := os.Open(gtf)
	defer gtfF.Close()
	if err != nil {
		log.Fatal(err)
	}

	Genes := make(map[string]*genodatastruct.Gene)
	allgene := geneset[0] == "all"
	scanner := bufio.NewScanner(gtfF)
	//var duration time.Duration
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#") { //skip comment lines
			continue
		}
		//parsing the gtf record
		//start := time.Now()
		fields := strings.Split(line, "\t")
		attributes := parseAttr(fields[len(fields)-1])
		//duration += time.Since(start)
		gene, geneid := attributes["gene_name"], attributes["gene_id"]

		//transform gtf record into corresponding data structure, Gene, Transcript
		if FeatureType(fields[2]) == Geneline {
			if allgene {
				Genes[geneid] = makeGene(fields, attributes)
			} else if _, has := hasgene(geneset, gene); has {
				Genes[geneid] = makeGene(fields, attributes)
			} else { //skip the gene
				continue
			}
		}

		if FeatureType(fields[2]) == Transcriptline {
			if allgene {
				Genes[geneid].Transcripts = append(Genes[geneid].Transcripts, makeTranscript(fields, attributes))
			} else if _, has := hasgene(geneset, gene); has {
				Genes[geneid].Transcripts = append(Genes[geneid].Transcripts, makeTranscript(fields, attributes))
			}
		}

		if FeatureType(fields[2]) == Exonline {
			idx := len(Genes[geneid].Transcripts) - 1
			start, _ := strconv.Atoi(fields[3])
			end, _ := strconv.Atoi(fields[4])
			if allgene {
				Genes[geneid].Transcripts[idx].Exons = append(Genes[geneid].Transcripts[idx].Exons, genodatastruct.Coor{start, end})
			} else if _, has := hasgene(geneset, gene); has {
				Genes[geneid].Transcripts[idx].Exons = append(Genes[geneid].Transcripts[idx].Exons, genodatastruct.Coor{start, end})
			}
		}
	}
	//println(duration)
	return Genes
}

//Initialize Gene struct from gene record in gtf
func makeGene(fields []string, attributes map[string]string) *genodatastruct.Gene {
	//gtf gene line has been splited and provided as arguments
	//var result Gene
	gene := attributes["gene_name"]
	start, _ := strconv.Atoi(fields[3])
	end, _ := strconv.Atoi(fields[4])
	return &genodatastruct.Gene{
		GeneName:   gene,
		Chromosome: genodatastruct.ChroSym(fields[0]),
		Coordinate: genodatastruct.Coor{start, end},
		Strand:     fields[6],
		Attributes: attributes,
	}
}

func makeTranscript(fields []string, attributes map[string]string) *genodatastruct.Transcript {
	transcript := attributes["transcript_name"]
	start, _ := strconv.Atoi(fields[3])
	end, _ := strconv.Atoi(fields[4])
	return &genodatastruct.Transcript{
		TranscriptName: transcript,
		Chromosome:     genodatastruct.ChroSym(fields[0]),
		Coordinate:     genodatastruct.Coor{start, end},
		Strand:         fields[6],
		Attributes:     attributes,
	}
}

//parseAttr splite attribute string by ";" and create key/val map
func parseAttr(attribute string) map[string]string {
	attriMap := make(map[string]string, 25)
	start, key := 0, ""
	for i, c := range attribute {
		if c == ' ' && start <= i {
			key = attribute[start:i]
			start = i + 1
		} else if c == ';' {
			val := attribute[start:i]
			attriMap[key] = strings.Trim(val, "\"")
			start = i + 2
		}
	}
	return attriMap
}

//check whether an item in a set
func hasitem(set []string, item string) (int, bool) {
	for i, rec := range set {
		if rec == item {
			return i, true
		}
	}
	return -1, false
}

//check whether a gene in a gene set, alias for general method of hasitem
//var hasgene hasItem
var hasgene = hasitem

//parsegtf provide a fast way to process a gtf line into data structure
//by iterate the strings only one time. Alternatively, split into fields
//and then splite attribute roughly equals to 2.5 rounds of iteration over the string
//func parsegtf(line string) (fields []string, attriMap map[string]string) {
//	//attriMap = make(map[string]string)
//	fields = make([]string, 8)
//	start, fieldCnt := 0, 0
//	for i, c := range line {
//		if c == '\t' {
//			fields[fieldCnt] = line[start:i]
//			start = i+1
//			fieldCnt++
//		}
//		if fieldCnt == 8 { //all left will be attribute
//			attriMap = parseAttrfast(line[start:])
//			break
//		}
//	}
//	return fields, attriMap
//}

//parseAttr splite attribute string by ";" and create key/val map
//func parseAttr(attribute string) map[string]string {
//	attriMap := make(map[string]string, 25)
//	fields := strings.Split(strings.TrimSuffix(attribute, ";"), ";")
//	//split each key/val pair of attri field
//	for _, pair := range fields {
//		kv := strings.Split(strings.TrimPrefix(pair, " "), " ")
//		if len(kv) != 2 {
//			log.Fatal("Ill-formated attribute field")
//		}
//		key, val := kv[0], strings.Trim(kv[1], "\"")
//		attriMap[key] = val
//	}
//	return attriMap
//}

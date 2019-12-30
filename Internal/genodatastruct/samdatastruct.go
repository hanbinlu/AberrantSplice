package genodatastruct

import (
	"log"
	"strconv"
	"strings"
)

//Only take sam fields related to mapping info
//not full support for all the sam fields yet
type SamRecPartial struct {
	Flag       int64
	CIGAR      string
	Pos        int
	Chromasome string
	MAPQ       int
}

func (sr *SamRecPartial) Strand() string {
	bitwise := strconv.FormatInt(sr.Flag, 2)
	if string(bitwise[len(bitwise)-5]) == "1" {
		return "-"
	} else {
		return "+"
	}
}

//Sample sam alignment record
//00953:17:HTJKFDSXX:2:1101:8287:1031    153     8       22663793        60      25M1671N76M     =       22663793        0
//CACTCACGGGACCTCGGACTTTGCTACAGGCGATCTTCAGGAGGTTCCAGAGTTCCTTCTGCCGCTTCTCCTGAAGCCGGACCACGGTCCTCTCGTCCTCG
//FFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF,F:FFFFF:
//AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:101        YT:Z:UP XS:A:-  NH:i:1

//CigarOp represent the following method. Return the corresponding regions
//of the reference (genome coordinate)
//+----+---------------------------------------------------------------+
//| op | Description                                                   |
//+----+---------------------------------------------------------------+
//| M  | Alignment match (can be a sequence match or mismatch)         |
//+----+---------------------------------------------------------------+
//| I  | Insertion to the reference                                    |
//+----+---------------------------------------------------------------+
//| D  | Deletion from the reference                                   |
//+----+---------------------------------------------------------------+
//| N  | Skipped region from the reference                             |
//+----+---------------------------------------------------------------+
//| S  | Soft clip on the read (clipped sequence present in <seq>)     |
//+----+---------------------------------------------------------------+
//| H  | Hard clip on the read (clipped sequence NOT present in <seq>) |
//+----+---------------------------------------------------------------+
//| P  | Padding (silent deletion from the padded reference sequence)  |
//+----+---------------------------------------------------------------+
//| =  | Sequence match                                                |
//+----+---------------------------------------------------------------+
//| X  | Sequence mismatch                                             |
//+----+---------------------------------------------------------------+
var cigarop = map[string]func(int, int) Coor{
	"M": func(pos, val int) Coor { return Coor{pos, pos + val - 1} },
	"I": func(pos, val int) Coor { return Coor{pos - 1, pos - 1} },
	"D": func(pos, val int) Coor { return Coor{pos, pos + val - 1} },
	"N": func(pos, val int) Coor { return Coor{pos, pos + val - 1} },
	"S": func(pos, val int) Coor { return Coor{pos - 1, pos - 1} },
	//"H": func(pos, val int) Coor { return Coor{pos, pos + val - 1} },
	//"=": func(pos, val int) Coor { return Coor{pos, pos + val - 1} },
	//"X": func(pos, val int) Coor { return Coor{pos, pos + val - 1} },
}

type op struct {
	pos   int
	steps int
	op    string
}

func (m *op) refreg() Coor {
	_, ok := cigarop[m.op]
	if !ok {
		log.Fatalln("Does not hard coded operator ", m.op)
	}
	return cigarop[m.op](m.pos, m.steps)
}

//RegionAligned parses CIGAR string and return the aligned region
//of the reference genome. Especially if intron exists, return
//more than one segment
func (sr SamRecPartial) RegionAligned() []Coor {
	//take each pattern of num/op in cigar
	//and walk on the reference instructed by num/op
	//to generate aligned region
	digits := "0123456789"
	start := 0
	aligned := []Coor{}
	walkfrom := sr.Pos
	for i, c := range sr.CIGAR {
		if !strings.ContainsRune(digits, c) { //trigger operation
			n, _ := strconv.Atoi(sr.CIGAR[start:i])
			operator := op{walkfrom, n, string(sr.CIGAR[i])}
			extension := operator.refreg()
			//combine extension to walked regions
			if operator.op == "N" { //intron contained splited read
				//initiate a new segment that skipping intron
				aligned = append(aligned, Coor{extension.End + 1, extension.End + 1})
			} else {
				lastSeg := &aligned[len(aligned)-1]
				//Continueous and updatable
				if extension.Start-lastSeg.End <= 1 && extension.End > lastSeg.End {
					lastSeg.End = extension.End
				} else {
					log.Fatalln("Broken segments")
				}
			}
			walkfrom = extension.End + 1
			start = i + 1
		} else {
			continue
		}
	}
	return aligned
}

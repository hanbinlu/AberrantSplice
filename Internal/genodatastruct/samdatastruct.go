package genodatastruct

import (
	"log"
	"strconv"
	"strings"
)

type SamFlag int

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

//CIGAR struct have data of mapping start position and cigar strings
//Methods will return the alignment scheme from the CIGAR
type CIGAR struct {
	cigar string
	pos   int
	chro  string
}

//RegionAligned parses CIGAR string and return the aligned region
//of the reference genome. Especially if intron exists, return
//more than one segment
func (cigar CIGAR) RegionAligned() []Coor {
	//take each pattern of num/op in cigar
	//and walk on the reference instructed by num/op
	//to generate aligned region
	digits := "0123456789"
	start := 0
	aligned := []Coor{}
	walkfrom := cigar.pos
	for i, c := range cigar.cigar {
		if !strings.ContainsRune(digits, c) { //trigger operation
			n, _ := strconv.Atoi(cigar.cigar[start:i])
			operator := op{walkfrom, n, string(cigar.cigar[i])}
			extension := operator.refreg()
			//combine extension to walked regions
			if operator.op == "N" { //intron contained splited read
				//add intron segment and a new segment
				aligned = append(aligned, extension, Coor{extension.End + 1, extension.End + 1})
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

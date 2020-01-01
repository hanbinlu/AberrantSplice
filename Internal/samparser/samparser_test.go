package samparser

import (
	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
	"os"
	"testing"
)

func TestParseSam(t *testing.T) {
	sam := os.Args[1]
	mapqfilter := func(s genodatastruct.SamRecPartial) bool {
		if s.MAPQ > 10 {
			return true
		} else {
			return false
		}
	}
	samchan := ParseSam(sam, mapqfilter)
	cnt := 0
	total := 0
	for i := range samchan {
		total++
		if i.Chromosome == "1" {
			cnt++
		}
	}
	println(total, "#####", cnt)
}

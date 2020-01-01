package samparser

import (
	"bufio"
	"github.com/Hanbin/AberrantSplice/Internal/genodatastruct"
	"log"
	"os"
	"strconv"
	"strings"
)

//ParseSam filter reads and parse them into SamRecPartial struct
//and generating a channel of iterator
func ParseSam(sam string, filter func(genodatastruct.SamRecPartial) bool) <-chan genodatastruct.SamRecPartial {
	out := make(chan genodatastruct.SamRecPartial, 100)
	go func() {
		samF, err := os.Open(sam)
		if err != nil {
			log.Fatal(err)
		}
		scanner := bufio.NewScanner(samF)
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, "@") {
				continue
			}
			fields := strings.Split(line, "\t")
			flag, _ := strconv.Atoi(fields[1])
			if flag == 4 { //unmapped
				continue
			}
			mapq, _ := strconv.Atoi(fields[4])
			pos, _ := strconv.Atoi(fields[3])
			temp := genodatastruct.SamRecPartial{
				Flag:       int64(flag),
				MAPQ:       mapq,
				Pos:        pos,
				Chromosome: genodatastruct.ChroSym(fields[2]),
				CIGAR:      fields[5],
			}
			if filter(temp) {
				out <- temp
			}
		}
		close(out)
	}()
	return out
}

package main

import (
	"flag"
	"fmt"
	"log"
)

var (
	inputFile  = flag.String("i", "", "Input file (required)")
	outputFile = flag.String("o", "", "Output file (required)")
	adapter    = flag.String("a", "", "Adapter sequence (required)")
	minLen     = flag.Int("minLen", 18, "Minimum length of read")
	trim5      = flag.Int("trim5", 0, "5' trim length")
	trim3      = flag.Int("trim3", 0, "3' trim length")
	min5Match  = flag.Int("min5Match", 8, "Minimum match length at 5' end")
	maxError   = flag.Float64("maxError", 0.1, "Maximum mean error rate")
)

func main() {
	flag.Parse()

	if *inputFile == "" || *outputFile == "" || *adapter == "" {
		fmt.Println("Missing required arguments")
		flag.Usage()
		return
	}

	err := ProcessReadsFast(*inputFile, *outputFile, *adapter, *minLen, *trim5, *trim3, *min5Match, *maxError)

	if err != nil {
		log.Fatalf("Error processing reads: %v", err)
	} else {
		fmt.Println("\nTrimming completed")
	}
}

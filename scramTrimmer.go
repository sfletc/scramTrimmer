package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"time"

	"github.com/fatih/color"
	"github.com/klauspost/pgzip"
)

type FastqRead struct {
	Header   string
	Sequence string
	Quality  string
}

// Rest of the utility functions remain the same
func phred33ToError(qual byte) float64 {
	return math.Pow(10, -(float64(qual)-33)/10.0)
}

func meanError(quality []byte) float64 {
	total := 0.0
	for _, q := range quality {
		total += phred33ToError(q)
	}
	return total / float64(len(quality))
}

func trimRead(read *FastqRead, adapter string, minLen, trim5, trim3, min5Match int, maxError float64) (*FastqRead, error) {
	adapterIndex := strings.Index(read.Sequence, adapter[:min5Match])

	if adapterIndex == -1 {
		return nil, fmt.Errorf("adapter missing")
	}

	start := trim5
	end := adapterIndex - trim3

	if end-start < minLen {
		return nil, fmt.Errorf("too short")
	}

	trimmedSequence := read.Sequence[start:end]
	trimmedQuality := read.Quality[start:end]

	if meanError([]byte(trimmedQuality)) >= maxError {
		return nil, fmt.Errorf("low quality")
	}

	trimmedRead := &FastqRead{
		Header:   read.Header,
		Sequence: trimmedSequence,
		Quality:  trimmedQuality,
	}
	return trimmedRead, nil
}

// Channel-based batch processor
func processBatch(
	batch []*FastqRead,
	adapter string,
	minLen, trim5, trim3, min5Match int,
	maxError float64,
	resultsChan chan<- *FastqRead,
	wg *sync.WaitGroup,
	adapterMissingCount, tooShortCount, lowQualityCount *int64,
) {
	defer wg.Done()

	for _, read := range batch {
		trimmedRead, err := trimRead(read, adapter, minLen, trim5, trim3, min5Match, maxError)
		if err != nil {
			switch err.Error() {
			case "adapter missing":
				atomic.AddInt64(adapterMissingCount, 1)
			case "too short":
				atomic.AddInt64(tooShortCount, 1)
			case "low quality":
				atomic.AddInt64(lowQualityCount, 1)
			}
			continue
		}
		resultsChan <- trimmedRead
	}
}

// Writer goroutine
func writeResults(
	writer *bufio.Writer,
	resultsChan <-chan *FastqRead,
	doneChan chan<- struct{},
	totalTrimmedReads *int64,
) {
	for read := range resultsChan {
		writer.WriteString(read.Header + "\n")
		writer.WriteString(read.Sequence + "\n")
		writer.WriteString("+\n")
		writer.WriteString(read.Quality + "\n")
		atomic.AddInt64(totalTrimmedReads, 1)
	}
	writer.Flush()
	close(doneChan)
}

func Comma(value int64) string {
	str := strconv.FormatInt(value, 10)
	result := ""
	count := 0
	for i := len(str) - 1; i >= 0; i-- {
		if count > 0 && count%3 == 0 {
			result = "," + result
		}
		result = string(str[i]) + result
		count++
	}
	return result
}

func ProcessReadsFast(inputFile, outputFile, adapter string, minLen, trim5, trim3, min5Match int, maxError float64) error {
	startTime := time.Now()

	inFile, err := os.Open(inputFile)
	if err != nil {
		return err
	}
	defer inFile.Close()

	gr, err := pgzip.NewReader(inFile)
	if err != nil {
		return err
	}
	defer gr.Close()

	outFile, err := os.Create(outputFile)
	if err != nil {
		return err
	}
	defer outFile.Close()

	gw := pgzip.NewWriter(outFile)
	defer gw.Close()
	writer := bufio.NewWriter(gw)

	// Create channels for processing
	resultsChan := make(chan *FastqRead, 1000) // Buffer size can be adjusted
	doneChan := make(chan struct{})

	var wg sync.WaitGroup
	var adapterMissingCount, tooShortCount, lowQualityCount int64
	var totalReads, totalTrimmedReads int64

	// Start writer goroutine
	go writeResults(writer, resultsChan, doneChan, &totalTrimmedReads)

	const batchSize = 10000 // Smaller batch size for better memory management
	scanner := bufio.NewScanner(gr)
	reads := make([]*FastqRead, 0, batchSize)

	// Process reads in batches
	for scanner.Scan() {
		header := scanner.Text()
		if !strings.HasPrefix(header, "@") {
			return fmt.Errorf("invalid fastq file: expected '@' at the beginning of header line, got: %s", header)
		}

		scanner.Scan()
		sequence := scanner.Text()

		scanner.Scan()
		plus := scanner.Text()
		if plus != "+" {
			return fmt.Errorf("invalid fastq file: expected '+' line, got: %s", plus)
		}

		scanner.Scan()
		quality := scanner.Text()
		if len(sequence) != len(quality) {
			return fmt.Errorf("invalid fastq file: sequence and quality strings must have the same length, got: %d and %d", len(sequence), len(quality))
		}

		reads = append(reads, &FastqRead{
			Header:   header,
			Sequence: sequence,
			Quality:  quality,
		})
		totalReads++

		if len(reads) == batchSize {
			wg.Add(1)
			go processBatch(reads, adapter, minLen, trim5, trim3, min5Match, maxError, resultsChan, &wg, &adapterMissingCount, &tooShortCount, &lowQualityCount)
			reads = make([]*FastqRead, 0, batchSize)
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error reading file: %v", err)
	}

	// Process remaining reads
	if len(reads) > 0 {
		wg.Add(1)
		go processBatch(reads, adapter, minLen, trim5, trim3, min5Match, maxError, resultsChan, &wg, &adapterMissingCount, &tooShortCount, &lowQualityCount)
	}

	// Wait for all processing to complete
	wg.Wait()
	close(resultsChan)

	// Wait for writer to finish
	<-doneChan

	// Calculate final statistics
	trimmedReadPercentage := (float64(totalTrimmedReads) / float64(totalReads)) * 100

	duration := time.Since(startTime)
	fmt.Printf("\nTotal reads: %s\n", Comma(totalReads))
	fmt.Printf("Trimmed reads: %s\n", Comma(totalTrimmedReads))
	color.HiGreen("Percentage of trimmed reads: %.2f%%\n", trimmedReadPercentage)
	color.HiMagenta("\nAdapter missing count: %s\n", Comma(adapterMissingCount))
	color.HiMagenta("Too short count: %s\n", Comma(tooShortCount))
	color.HiMagenta("Low quality count: %s\n", Comma(lowQualityCount))
	fmt.Printf("\nApplication execution time: %s\n", duration)

	return nil
}

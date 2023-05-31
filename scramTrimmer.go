package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io/ioutil"
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

type Adapter struct {
	Seq string
}

/*
phred33ToError is a function that converts a single Phred quality score encoded as a Phred+33 ASCII value into an error probability.

Parameter:
- qual: A byte representing a Phred+33 ASCII-encoded quality score. In a FASTQ file,
  quality scores are encoded as ASCII characters with Phred+33 encoding (also known as Sanger encoding),
  where the character's ASCII value minus 33 gives the Phred quality score.

The function subtracts 33 from the ASCII value of the input byte to derive the Phred quality score.
This score represents the base-calling error probabilities on a logarithmic scale.
This function then converts this logarithmic score into an error probability using the formula:

  error = 10 ^ (-(Phred score) / 10)

which is the standard formula to convert Phred scores to error probabilities.

Returns:
- A float64 representing the error probability associated with the provided Phred+33 ASCII-encoded quality score.
*/

func phred33ToError(qual byte) float64 {
	return math.Pow(10, -(float64(qual)-33)/10.0)
}

/*
meanError is a function that calculates the mean error probability of a sequence read based on its quality scores.

Parameter:
- quality: A byte slice representing Phred+33 ASCII-encoded quality scores. Each byte in this slice corresponds
  to a base in the sequence read, with the byte's value representing the error probability of that base call
  being incorrect.

The function iterates over the byte slice and converts each Phred+33 ASCII-encoded quality score into an error
probability using the phred33ToError function. It then sums up these error probabilities and divides the sum
by the total number of quality scores to calculate the mean error probability.

Returns:
- A float64 representing the mean error probability of the sequence read.
*/

func meanError(quality []byte) float64 {
	total := 0.0
	for _, q := range quality {
		total += phred33ToError(q)
	}
	return total / float64(len(quality))
}

/*
trimRead is a function that trims a FASTQ sequence read based on the provided parameters.

Parameters:
- read: A pointer to a FastqRead object. This object represents a single read from a FASTQ file,
  which includes a header, sequence and a quality score for each base in the sequence.
- adapter: A string representing the sequence of the adapter that is to be removed from the read.
- minLen: An integer specifying the minimum length a read can be after trimming.
- trim5: An integer specifying the number of bases to be removed from the 5' end of the read.
- trim3: An integer specifying the number of bases to be removed from the 3' end of the read after adapter trimming.
- min5Match: An integer specifying the minimum number of bases of the 5' end of
  the adapter that match the read to be considered present.
- maxError: A float64 value that represents the maximum average base calling error rate allowed for
  the trimmed read.

The function first checks whether the adapter is present in the read. If it is not, an error is returned.

The function then checks if the length of the trimmed read will be less than the specified minimum length.
If it is, an error is returned.

The read is trimmed according to the provided trim5, trim3 and adapter parameters.

The quality of the read is then assessed by calculating the mean base calling error rate.
If this error rate is greater than or equal to the maximum allowed error rate, an error is returned.

If all checks pass, a new FastqRead object is created with the header, trimmed sequence and quality score,
and it is returned along with nil error.

Returns:
- A pointer to a FastqRead object containing the header, trimmed sequence and its corresponding quality score.
  If any of the checks fail, nil is returned instead.
- An error that is nil if all checks pass. If any of the checks fail, an appropriate error message is returned.
*/

func trimRead(read *FastqRead, adapter string, minLen, trim5, trim3, min5Match int, maxError float64) (*FastqRead, error) {
	adapterIndex := strings.Index(read.Sequence, adapter[:min5Match])

	// Check if adapter is present with min5Match
	if adapterIndex == -1 {
		return nil, fmt.Errorf("adapter missing")
	}

	start := trim5
	end := adapterIndex - trim3

	// Check if read is long enough after trimming
	if end-start < minLen {
		return nil, fmt.Errorf("too short")
	}

	// Trim the read
	trimmedSequence := read.Sequence[start:end]
	trimmedQuality := read.Quality[start:end]

	// Check the mean error
	if meanError([]byte(trimmedQuality)) >= maxError {
		return nil, fmt.Errorf("low quality")
	}

	// Return the trimmed read if all checks pass
	trimmedRead := &FastqRead{
		Header:   read.Header,
		Sequence: trimmedSequence,
		Quality:  trimmedQuality,
	}
	return trimmedRead, nil
}

const batchSize = 1_000_000

// Read the entire input file into memory
func readFile(inputFile string) ([]byte, error) {
	return ioutil.ReadFile(inputFile)
}

// Prepare to write the output file
func createOutputFile(outputFile string) (*os.File, error) {
	return os.Create(outputFile)
}

/*
processBatch is a function that trims a batch of FASTQ sequence reads concurrently, and updates the
statistics on reads with missing adapters, those that are too short after trimming, and those with low quality.

Parameters:
- batch: A slice of pointers to FastqRead objects, each representing a single read from a FASTQ file.
- adapter: A string representing the sequence of the adapter that is to be removed from the reads.
- minLen: An integer specifying the minimum length a read can be after trimming.
- trim5: An integer specifying the number of bases to be removed from the 5' end of the read.
- trim3: An integer specifying the number of bases to be removed from the 3' end of the read after adapter trimming.
- min5Match: An integer specifying the minimum number of bases of the 5' end of
  the adapter that match the read to be considered present.
- maxError: A float64 value that represents the maximum average base calling error rate allowed for
  the trimmed read.
- wg: A pointer to a sync.WaitGroup object to allow the parent goroutine to wait for the completion of this function.
- mutex: A pointer to a sync.Mutex object to ensure safe concurrent access to shared data, specifically the slice of trimmed reads.
- trimmedReads: A pointer to a slice of pointers to FastqRead objects to hold the trimmed reads.
- adapterMissingCount: A pointer to an int64 to keep count of the number of reads with missing adapters.
- tooShortCount: A pointer to an int64 to keep count of the number of reads that are too short after trimming.
- lowQualityCount: A pointer to an int64 to keep count of the number of reads with low quality.

The function loops over the batch of reads, trims each read using the trimRead function, and checks for any errors.
If an error occurs, the function identifies the type of error and updates the corresponding error count using
atomic operations to ensure thread safety.

If no error occurs, the trimmed read is safely added to the slice of trimmed reads using the provided mutex
to ensure thread safety.

The function signals the completion of its execution to the wait group using the defer statement, ensuring
that it will be executed even if an error occurs during processing.

This function is designed to be run as a goroutine for concurrent processing of large FASTQ files.
*/

func processBatch(batch []*FastqRead, adapter string, minLen, trim5, trim3, min5Match int, maxError float64, wg *sync.WaitGroup, mutex *sync.Mutex, trimmedReads *[]*FastqRead, adapterMissingCount, tooShortCount, lowQualityCount *int64) {
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

		mutex.Lock()
		*trimmedReads = append(*trimmedReads, trimmedRead)
		mutex.Unlock()
	}
}

/*
writeTrimmedReads is a function that writes a list of trimmed sequence reads in FASTQ format to a buffered writer.

Parameters:
- trimmedReads: A slice of pointers to FastqRead objects. Each object represents a single trimmed sequence read,
  including the read's header, sequence, and quality score.
- writer: A pointer to a bufio.Writer object. This writer is used to write the trimmed reads to an output stream
  (like a file) in a buffered manner, which can improve performance when writing large amounts of data.

The function iterates over the list of trimmed reads and writes the components of each read to the writer
in the standard FASTQ format, which is:

@HEADER
SEQUENCE
+
QUALITY

The function handles errors that may occur during writing. If an error occurs while writing any component
of a read, the function immediately returns the error.

Returns:
- An error which will be nil if all writes are successful. If a write fails, the function returns the error
  associated with that write.
*/

func writeTrimmedReads(trimmedReads []*FastqRead, writer *bufio.Writer) error {
	for _, read := range trimmedReads {
		_, err := writer.WriteString(read.Header + "\n")
		if err != nil {
			return err
		}
		_, err = writer.WriteString(read.Sequence + "\n")
		if err != nil {
			return err
		}
		_, err = writer.WriteString("+\n")
		if err != nil {
			return err
		}
		_, err = writer.WriteString(read.Quality + "\n")
		if err != nil {
			return err
		}
	}
	return nil
}

/*
Comma is a function that formats an int64 value as a string with commas as thousand separators.

Parameter:
- value: An int64 value that needs to be converted into a string and formatted with commas.

The function first converts the int64 value to a string. It then traverses the string in reverse
(from right to left) and appends each digit to the result string. For every three digits, it also
prepends a comma to the result string. This results in a string representation of the number
with commas as thousand separators.

Returns:
- A string representing the input value formatted with commas as thousand separators.
*/

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

/*
ProcessReadsFast is a function that processes sequence reads from a compressed FASTQ file in a parallel manner, trims adapters, checks the quality, and writes the trimmed reads to a new compressed FASTQ file. It also reports statistics on the reads processed.

Parameters:
- inputFile: A string representing the path to the input FASTQ file.
- outputFile: A string representing the path to the output FASTQ file.
- adapter: A string representing the sequence of the adapter that needs to be removed from the reads.
- minLen: An integer specifying the minimum length a read can be after trimming.
- trim5: An integer specifying the number of bases to be removed from the 5' end of the read.
- trim3: An integer specifying the number of bases to be removed from the 3' end of the read after adapter trimming.
- min5Match: An integer specifying the minimum number of bases of the 5' end of
  the adapter that match the read to be considered present.
- maxError: A float64 value that represents the maximum average base calling error rate allowed for
  the trimmed read.

The function reads the input file, decompresses it using the pgzip package, and then scans it line by line to extract the FASTQ reads. Reads are processed in batches using the processBatch function, which is run as a goroutine for concurrency.

Counts are kept for total reads, trimmed reads, and the number of reads missing the adapter, that are too short after trimming, and that have low quality.

Upon completion of all processing, the function writes the trimmed reads to the output file using the writeTrimmedReads function and reports the processing statistics.

Returns:
- An error which will be nil if all operations are successful. If an error occurs, the function returns the error associated with that operation.
*/

func ProcessReadsFast(inputFile, outputFile, adapter string, minLen, trim5, trim3, min5Match int, maxError float64) error {
	startTime := time.Now()

	inFileBytes, err := readFile(inputFile)
	if err != nil {
		return err
	}
	inFile := bytes.NewReader(inFileBytes)

	gr, err := pgzip.NewReader(inFile)
	if err != nil {
		return err
	}
	defer gr.Close()

	outFile, err := createOutputFile(outputFile)
	if err != nil {
		return err
	}
	defer outFile.Close()

	gw := pgzip.NewWriter(outFile)
	defer gw.Close()
	writer := bufio.NewWriter(gw)

	scanner := bufio.NewScanner(gr)
	var wg sync.WaitGroup
	var adapterMissingCount, tooShortCount, lowQualityCount int64
	var totalReads, totalTrimmedReads int64 // New counters for total reads and trimmed reads
	var mutex sync.Mutex
	var trimmedReads []*FastqRead
	reads := make([]*FastqRead, 0, batchSize)

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
		totalReads++ // Increment total reads counter

		if len(reads) == batchSize {
			wg.Add(1)
			go processBatch(reads, adapter, minLen, trim5, trim3, min5Match, maxError, &wg, &mutex, &trimmedReads, &adapterMissingCount, &tooShortCount, &lowQualityCount)
			reads = make([]*FastqRead, 0, batchSize)
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error reading file: %v", err)
	}

	// process any remaining reads
	if len(reads) > 0 {
		wg.Add(1)
		go processBatch(reads, adapter, minLen, trim5, trim3, min5Match, maxError, &wg, &mutex, &trimmedReads, &adapterMissingCount, &tooShortCount, &lowQualityCount)
	}

	wg.Wait()

	totalTrimmedReads = int64(len(trimmedReads))                                      // Set total trimmed reads counter
	trimmedReadPercentage := (float64(totalTrimmedReads) / float64(totalReads)) * 100 // Calculate percentage of trimmed reads

	err = writeTrimmedReads(trimmedReads, writer)
	if err != nil {
		return err
	}
	writer.Flush()
	gw.Close()

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

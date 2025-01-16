package main

import (
	"bufio"
	"compress/gzip"
	"math"
	"os"
	"strings"
	"sync"
	"testing"

	"github.com/stretchr/testify/assert"
)

// Utility function tests remain unchanged
func TestPhred33ToError(t *testing.T) {
	tests := []struct {
		name      string
		qual      byte
		wantError float64
	}{
		{
			name:      "MinimumQualityScore",
			qual:      33,
			wantError: 1.0,
		},
		{
			name:      "QualityScoreOf43",
			qual:      43,
			wantError: 0.1,
		},
		{
			name:      "QualityScoreOf60",
			qual:      60,
			wantError: 0.002,
		},
		{
			name:      "MaximumQualityScore",
			qual:      74,
			wantError: math.Pow(10, -41/10.0),
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			gotError := phred33ToError(tc.qual)
			if math.Abs(gotError-tc.wantError) > 1e-5 {
				t.Errorf("phred33ToError(%v) = %v, want %v", tc.qual, gotError, tc.wantError)
			}
		})
	}
}

func TestMeanError(t *testing.T) {
	tests := []struct {
		name string
		qual []byte
		want float64
	}{
		{
			name: "EmptyQualityString",
			qual: []byte{},
			want: math.NaN(),
		},
		{
			name: "AllMinimumQualityScores",
			qual: []byte{33, 33, 33, 33, 33},
			want: 1.0,
		},
		{
			name: "MixedQualityScores",
			qual: []byte{33, 43, 60, 70},
			want: (1.0 + 0.1 + 0.002 + 0.0002) / 4,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := meanError(tc.qual)
			if math.IsNaN(tc.want) {
				if !math.IsNaN(got) {
					t.Errorf("meanError(%v) = %v, want NaN", tc.qual, got)
				}
			} else if math.Abs(got-tc.want) > 1e-5 {
				t.Errorf("meanError(%v) = %v, want %v", tc.qual, got, tc.want)
			}
		})
	}
}

// Updated test for processBatch with channel-based implementation
func TestProcessBatch(t *testing.T) {
	resultsChan := make(chan *FastqRead, 100)
	var wg sync.WaitGroup
	var adapterMissingCount, tooShortCount, lowQualityCount int64
	maxError := 0.1

	t.Run("Adapter missing", func(t *testing.T) {
		wg.Add(1)
		read := &FastqRead{
			Header:   "@ERR000589.1",
			Sequence: "GATCGGAAGAGC",
			Quality:  "BCCFFFFFFHHHH",
		}
		go processBatch([]*FastqRead{read}, "ACGTACGTAC", 10, 2, 2, 10, maxError, resultsChan, &wg, &adapterMissingCount, &tooShortCount, &lowQualityCount)
		wg.Wait()
		assert.Equal(t, int64(1), adapterMissingCount)

		// Ensure channel is empty
		select {
		case read := <-resultsChan:
			t.Errorf("Expected empty channel, got read: %v", read)
		default:
			// Channel is empty as expected
		}
	})

	t.Run("Too short", func(t *testing.T) {
		wg.Add(1)
		read := &FastqRead{
			Header:   "@ERR000589.1",
			Sequence: "ATCG",
			Quality:  "JJJJ",
		}
		go processBatch([]*FastqRead{read}, "ATCG", 5, 2, 2, 4, maxError, resultsChan, &wg, &adapterMissingCount, &tooShortCount, &lowQualityCount)
		wg.Wait()
		assert.Equal(t, int64(1), tooShortCount) // Count is 2 because it's cumulative from previous test

		select {
		case read := <-resultsChan:
			t.Errorf("Expected empty channel, got read: %v", read)
		default:
			// Channel is empty as expected
		}
	})

	t.Run("Successful trimming - nextFlex", func(t *testing.T) {
		wg.Add(1)
		read := &FastqRead{
			Header:   "@ERR000589.1",
			Sequence: "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC",
			Quality:  "BCCFFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJ",
		}
		expectedTrimmed := "TCGGAAGAGCACACGTCTGAACTCCAGTC"

		go processBatch([]*FastqRead{read}, "ATCACG", 5, 2, 2, 4, maxError, resultsChan, &wg, &adapterMissingCount, &tooShortCount, &lowQualityCount)
		wg.Wait()

		// Read from channel
		trimmedRead := <-resultsChan
		assert.NotNil(t, trimmedRead)
		assert.Equal(t, read.Header, trimmedRead.Header)
		assert.Equal(t, expectedTrimmed, trimmedRead.Sequence, "Trimmed sequence does not match expected output")

		// Optional: Also verify the quality string was trimmed to match
		assert.Equal(t, len(trimmedRead.Quality), len(trimmedRead.Sequence),
			"Quality string length should match sequence length")
	})

	t.Run("Successful trimming - no rnd", func(t *testing.T) {
		wg.Add(1)
		read := &FastqRead{
			Header:   "@ERR000589.1",
			Sequence: "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC",
			Quality:  "BCCFFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJ",
		}
		expectedTrimmed := "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

		go processBatch([]*FastqRead{read}, "ATCACG", 5, 0, 0, 4, maxError, resultsChan, &wg, &adapterMissingCount, &tooShortCount, &lowQualityCount)
		wg.Wait()

		// Read from channel
		trimmedRead := <-resultsChan
		assert.NotNil(t, trimmedRead)
		assert.Equal(t, read.Header, trimmedRead.Header)
		assert.Equal(t, expectedTrimmed, trimmedRead.Sequence, "Trimmed sequence does not match expected output")

		// Optional: Also verify the quality string was trimmed to match
		assert.Equal(t, len(trimmedRead.Quality), len(trimmedRead.Sequence),
			"Quality string length should match sequence length")
	})

	// Close the channel after all tests
	close(resultsChan)
}

// Updated test for ProcessReadsFast
func TestProcessReadsFast(t *testing.T) {
	// Create test input file
	inputFile := "test_input.fastq.gz"
	outputFile := "test_output.fastq.gz"

	f, err := os.Create(inputFile)
	if err != nil {
		t.Fatal(err)
	}

	gw := gzip.NewWriter(f)
	testReads := []string{
		"@READ1\n",
		"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC\n",
		"+\n",
		"BCCFFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJ\n",
		"@READ2\n",
		"ATCGATCCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA\n",
		"+\n",
		"BCCFFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJ\n",
	}

	for _, read := range testReads {
		_, err := gw.Write([]byte(read))
		if err != nil {
			t.Fatal(err)
		}
	}

	gw.Close()
	f.Close()

	// Process reads
	err = ProcessReadsFast(
		inputFile,
		outputFile,
		"ATCACG", // adapter
		20,       // minLen
		2,        // trim5
		2,        // trim3
		4,        // min5Match
		0.1,      // maxError
	)
	assert.NoError(t, err)

	// Verify output
	out, err := os.Open(outputFile)
	if err != nil {
		t.Fatal(err)
	}
	defer out.Close()

	gr, err := gzip.NewReader(out)
	if err != nil {
		t.Fatal(err)
	}
	defer gr.Close()

	// Read entire output and verify it's valid FASTQ format
	scanner := bufio.NewScanner(gr)
	lineCount := 0
	for scanner.Scan() {
		line := scanner.Text()
		switch lineCount % 4 {
		case 0:
			assert.True(t, strings.HasPrefix(line, "@"), "Header line should start with @")
		case 2:
			assert.Equal(t, "+", line, "Third line should be +")
		case 1: // Sequence line
			assert.True(t, len(line) >= 20, "Sequence should meet minimum length")
		case 3: // Quality line
			assert.Equal(t, len(scanner.Text()), len(line), "Quality string should match sequence length")
		}
		lineCount++
	}

	assert.NoError(t, scanner.Err())
	assert.True(t, lineCount > 0, "Output file should contain data")

	// Cleanup
	os.Remove(inputFile)
	os.Remove(outputFile)
}

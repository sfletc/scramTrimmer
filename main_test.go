package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"reflect"
	"sync"
	"testing"

	"github.com/stretchr/testify/assert"
)

// TestPhred33ToError verifies that the function correctly calculates error probabilities from PHRED33 scores.
func TestPhred33ToError(t *testing.T) {
	// Define test cases
	tests := []struct {
		name      string
		qual      byte
		wantError float64
	}{
		{
			name:      "MinimumQualityScore",
			qual:      33, // Minimum possible PHRED33 score, which should result in an error probability of 1.0
			wantError: 1.0,
		},
		{
			name:      "QualityScoreOf43",
			qual:      43, // PHRED33 score of 40, which should result in an error probability of 0.0001
			wantError: 0.1,
		},
		{
			name:      "QualityScoreOf60",
			qual:      60, // PHRED33 score of 50, which should result in an error probability of 0.00001
			wantError: 0.002,
		},
		{
			name:      "MaximumQualityScore",
			qual:      74, // Maximum possible PHRED33 score (assuming the use of the extended PHRED33 scale), which should result in the smallest error probability
			wantError: math.Pow(10, -41/10.0),
		},
	}

	// Run test cases
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
			qual: []byte{}, // An empty quality string should return NaN, because division by zero.
			want: math.NaN(),
		},
		{
			name: "AllMinimumQualityScores",
			qual: []byte{33, 33, 33, 33, 33}, // All minimum quality scores should result in a mean error rate of 1.
			want: 1.0,
		},
		{
			name: "MixedQualityScores",
			qual: []byte{33, 43, 60, 70}, // A mix of quality scores should result in a calculated mean error rate.
			want: (1.0 + 0.1 + 0.002 + 0.0002) / 4,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := meanError(tc.qual)

			// Need to handle the special case where we expect NaN, since NaN != NaN.
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

// TestTrimRead verifies the trimRead function with various input scenarios
func TestTrimRead(t *testing.T) {
	tests := []struct {
		name         string
		read         *FastqRead
		adapter      string
		minLen       int
		trim5        int
		trim3        int
		min5Match    int
		maxError     float64
		expectedRead *FastqRead
		expectedErr  error
	}{
		{
			name: "AdapterMissing",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCG",
				Quality:  "FFFFFFFF",
			},
			adapter:      "TTTT",
			minLen:       4,
			trim5:        0,
			trim3:        0,
			min5Match:    3,
			maxError:     0.1,
			expectedRead: nil,
			expectedErr:  fmt.Errorf("adapter missing"),
		},
		{
			name: "ReadTooShort",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "FFFFFFFF",
			},
			adapter:      "ATCC",
			minLen:       10,
			trim5:        0,
			trim3:        0,
			min5Match:    4,
			maxError:     0.1,
			expectedRead: nil,
			expectedErr:  fmt.Errorf("too short"),
		},
		{
			name: "LowQuality",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "!!!!!!!!",
			},
			adapter:      "ATCC",
			minLen:       2,
			trim5:        0,
			trim3:        0,
			min5Match:    4,
			maxError:     0.1,
			expectedRead: nil,
			expectedErr:  fmt.Errorf("low quality"),
		},
		{
			name: "SuccessfulTrim",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "FFFFFFFF",
			},
			adapter:   "ATCC",
			minLen:    2,
			trim5:     0,
			trim3:     0,
			min5Match: 4,
			maxError:  0.1,
			expectedRead: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCG",
				Quality:  "FFFF",
			},
			expectedErr: nil,
		},
		{
			name: "SuccessfulTrimWithTrim5",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "FFFFFFFF",
			},
			adapter:   "ATCC",
			minLen:    2,
			trim5:     2,
			trim3:     0,
			min5Match: 4,
			maxError:  0.1,
			expectedRead: &FastqRead{
				Header:   "@Header",
				Sequence: "CG",
				Quality:  "FF",
			},
			expectedErr: nil,
		},
		{
			name: "SuccessfulTrimWithTrim3",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "FFFFFFFF",
			},
			adapter:   "ATCC",
			minLen:    2,
			trim5:     0,
			trim3:     2,
			min5Match: 4,
			maxError:  0.1,
			expectedRead: &FastqRead{
				Header:   "@Header",
				Sequence: "AT",
				Quality:  "FF",
			},
			expectedErr: nil,
		},
		{
			name: "SuccessfulTrimWithTrim5AndTrim3",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "FFFFFFFF",
			},
			adapter:   "ATCC",
			minLen:    2,
			trim5:     1,
			trim3:     1,
			min5Match: 4,
			maxError:  0.1,
			expectedRead: &FastqRead{
				Header:   "@Header",
				Sequence: "TC",
				Quality:  "FF",
			},
			expectedErr: nil,
		},
		{
			name: "Trim5LongerThanReadAfterAdapterTrimmingWithNewParameters",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "FFFFFFFF",
			},
			adapter:      "ATCC",
			minLen:       2,
			trim5:        5,
			trim3:        0,
			min5Match:    4,
			maxError:     0.1,
			expectedRead: nil,
			expectedErr:  fmt.Errorf("too short"),
		},
		{
			name: "Trim3LongerThanReadAfterAdapterTrimmingWithNewParameters",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "FFFFFFFF",
			},
			adapter:      "ATCC",
			minLen:       2,
			trim5:        0,
			trim3:        5,
			min5Match:    4,
			maxError:     0.1,
			expectedRead: nil,
			expectedErr:  fmt.Errorf("too short"),
		},
		{
			name: "Trim5AndTrim3LongerThanReadAfterAdapterTrimmingWithNewParameters",
			read: &FastqRead{
				Header:   "@Header",
				Sequence: "ATCGATCC",
				Quality:  "FFFFFFFF",
			},
			adapter:      "ATCC",
			minLen:       2,
			trim5:        3,
			trim3:        3,
			min5Match:    4,
			maxError:     0.1,
			expectedRead: nil,
			expectedErr:  fmt.Errorf("too short"),
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			gotRead, gotErr := trimRead(tc.read, tc.adapter, tc.minLen, tc.trim5, tc.trim3, tc.min5Match, tc.maxError)
			if !reflect.DeepEqual(gotRead, tc.expectedRead) {
				t.Errorf("got read %v, want read %v", gotRead, tc.expectedRead)
			}
			if gotErr != nil && tc.expectedErr != nil {
				if gotErr.Error() != tc.expectedErr.Error() {
					t.Errorf("got error %v, want error %v", gotErr, tc.expectedErr)
					return
				}
			} else if gotErr != tc.expectedErr {
				t.Errorf("got error %v, want error %v", gotErr, tc.expectedErr)
			}
		})
	}
}

func TestReadFile(t *testing.T) {
	// Create a temporary file
	tempFile, err := ioutil.TempFile("", "testfile")
	if err != nil {
		t.Fatalf("Failed to create temporary file: %v", err)
	}
	defer os.Remove(tempFile.Name()) // clean up

	testContent := "Test content for readFile function"
	_, err = tempFile.WriteString(testContent)
	if err != nil {
		t.Fatalf("Failed to write to temporary file: %v", err)
	}

	// Close the file (important to flush to disk)
	tempFile.Close()

	// Test readFile
	content, err := readFile(tempFile.Name())
	if err != nil {
		t.Fatalf("Failed to read file: %v", err)
	}

	// Assert the content is correct
	if string(content) != testContent {
		t.Fatalf("Expected file content to be '%s', got '%s'", testContent, string(content))
	}
}

func TestCreateOutputFile(t *testing.T) {
	// Temporary file for testing
	tempFile, err := ioutil.TempFile("", "testfile")
	if err != nil {
		t.Fatalf("Failed to create temporary file: %v", err)
	}
	defer os.Remove(tempFile.Name()) // clean up

	// Test createOutputFile
	file, err := createOutputFile(tempFile.Name())
	if err != nil {
		t.Fatalf("Failed to create output file: %v", err)
	}

	testContent := "Test content for createOutputFile function"
	_, err = file.WriteString(testContent)
	if err != nil {
		t.Fatalf("Failed to write to output file: %v", err)
	}

	// Close the file (important to flush to disk)
	file.Close()

	// Verify if content is written correctly
	content, err := ioutil.ReadFile(tempFile.Name())
	if err != nil {
		t.Fatalf("Failed to read file: %v", err)
	}

	// Assert the content is correct
	if string(content) != testContent {
		t.Fatalf("Expected file content to be '%s', got '%s'", testContent, string(content))
	}
}

func TestCreateOutputFileNoPermission(t *testing.T) {
	// Test creating a file in a directory without write permissions
	_, err := createOutputFile("/root/testfile")
	if err == nil {
		t.Fatalf("Expected error due to lack of write permissions, got nil")
	}
}

func TestProcessBatch2(t *testing.T) {
	mutex := &sync.Mutex{}
	trimmedReads := make([]*FastqRead, 0)
	adapterMissingCount, tooShortCount, lowQualityCount := int64(0), int64(0), int64(0)
	maxError := 0.1

	t.Run("Adapter missing", func(t *testing.T) {
		wg := &sync.WaitGroup{}
		wg.Add(1)
		read := &FastqRead{
			Header:   "@ERR000589.1 HSQ1004:134:C0D8DACXX:2:1101:1243:2213/1",
			Sequence: "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC",
			Quality:  "BCCFFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJ",
		}
		go processBatch([]*FastqRead{read}, "ACGTACGTAC", 10, 2, 2, 10, maxError, wg, mutex, &trimmedReads, &adapterMissingCount, &tooShortCount, &lowQualityCount)
		wg.Wait()
		assert.Equal(t, int64(1), adapterMissingCount)
	})

	t.Run("Too short", func(t *testing.T) {
		wg := &sync.WaitGroup{}
		wg.Add(1)
		read := &FastqRead{
			Header:   "@ERR000589.1 HSQ1004:134:C0D8DACXX:2:1101:1243:2213/1",
			Sequence: "ATCG",
			Quality:  "JJJJ",
		}
		go processBatch([]*FastqRead{read}, "ATCG", 5, 2, 2, 4, maxError, wg, mutex, &trimmedReads, &adapterMissingCount, &tooShortCount, &lowQualityCount)
		wg.Wait()
		assert.Equal(t, int64(1), tooShortCount)
	})

	t.Run("Low quality", func(t *testing.T) {
		wg := &sync.WaitGroup{}
		wg.Add(1)
		read := &FastqRead{
			Header:   "@ERR000589.1 HSQ1004:134:C0D8DACXX:2:1101:1243:2213/1",
			Sequence: "ATCGATCC",
			Quality:  "!!!!!!!!",
		}
		go processBatch([]*FastqRead{read}, "ATCC", 3, 0, 0, 4, maxError, wg, mutex, &trimmedReads, &adapterMissingCount, &tooShortCount, &lowQualityCount)
		wg.Wait()
		assert.Equal(t, int64(1), lowQualityCount)
	})

	t.Run("Successful trimming", func(t *testing.T) {
		wg := &sync.WaitGroup{}
		wg.Add(1)
		read := &FastqRead{
			Header:   "@ERR000589.1 HSQ1004:134:C0D8DACXX:2:1101:1243:2213/1",
			Sequence: "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC",
			Quality:  "BCCFFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJ",
		}
		go processBatch([]*FastqRead{read}, "ATCACG", 5, 2, 2, 4, maxError, wg, mutex, &trimmedReads, &adapterMissingCount, &tooShortCount, &lowQualityCount)
		wg.Wait()
		assert.Equal(t, 1, len(trimmedReads))
	})
}

func TestWriteTrimmedReads(t *testing.T) {
	reads := []*FastqRead{
		{
			Header:   "@SEQ_ID",
			Sequence: "ACTG",
			Quality:  "!!!!",
		},
		{
			Header:   "@SEQ_ID2",
			Sequence: "TGCA",
			Quality:  "****",
		},
	}

	buf := &bytes.Buffer{}
	writer := bufio.NewWriter(buf)

	err := writeTrimmedReads(reads, writer)
	assert.NoError(t, err)
	writer.Flush() // Ensure all data is written to the buffer

	expected := "@SEQ_ID\nACTG\n+\n!!!!\n@SEQ_ID2\nTGCA\n+\n****\n"
	assert.Equal(t, expected, buf.String())
}

func TestComma(t *testing.T) {
	tests := []struct {
		name     string
		input    int64
		expected string
	}{
		{
			name:     "Test 1 - Less than 1000",
			input:    123,
			expected: "123",
		},
		{
			name:     "Test 2 - Thousand",
			input:    1234,
			expected: "1,234",
		},
		{
			name:     "Test 3 - Million",
			input:    1234567,
			expected: "1,234,567",
		},
		{
			name:     "Test 4 - Billion",
			input:    1234567890,
			expected: "1,234,567,890",
		},
		{
			name:     "Test 5 - Trillion",
			input:    1234567890123,
			expected: "1,234,567,890,123",
		},
		{
			name:     "Test 6 - Zero",
			input:    0,
			expected: "0",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := Comma(tt.input)
			assert.Equal(t, tt.expected, result)
		})
	}
}

func TestProcessReadsFast(t *testing.T) {
	// Setup test data
	adapter := "ATCC"
	minLen := 4
	trim5 := 0
	trim3 := 0
	min5Match := 4
	maxError := 0.1

	// Setup test file paths
	inputFile := "test_input.fastq.gz"
	outputFile := "test_output.fastq.gz"

	// Create test input file
	f, err := os.Create(inputFile)
	if err != nil {
		t.Fatal(err)
	}

	gw := gzip.NewWriter(f)

	gw.Write([]byte("@READ1\n"))
	gw.Write([]byte("ATCGATCC\n"))
	gw.Write([]byte("+\n"))
	gw.Write([]byte("IIIIIIII\n"))

	gw.Close()
	f.Close()

	// Process reads
	err = ProcessReadsFast(inputFile, outputFile, adapter, minLen, trim5, trim3, min5Match, maxError)
	assert.NoError(t, err)

	// Open output file for validation
	out, err := os.Open(outputFile)
	if err != nil {
		t.Fatal(err)
	}
	defer out.Close()

	gr, err := gzip.NewReader(out)
	if err != nil {
		t.Fatal(err)
	}

	// Use bufio or ioutil package to read and validate the output
	scanner := bufio.NewScanner(gr)

	// Check first read
	assert.True(t, scanner.Scan())
	assert.Equal(t, "@READ1", scanner.Text())

	assert.True(t, scanner.Scan())
	// The trimmed sequence should be "ATCG" if the adapter was "ATCC" and we trimmed from the start of the adapter
	assert.Equal(t, "ATCG", scanner.Text())

	assert.True(t, scanner.Scan())
	assert.Equal(t, "+", scanner.Text())

	assert.True(t, scanner.Scan())
	// The quality string should match the length of the trimmed sequence
	assert.Equal(t, "IIII", scanner.Text())

	assert.False(t, scanner.Scan(), "Should be no more lines in the file")

	if scanner.Err() != nil {
		t.Errorf("Error scanning output file: %v", scanner.Err())
	}

	// Clean up test files
	os.Remove(inputFile)
	os.Remove(outputFile)
}

# RNA Adapter Trimmer

RNA Adapter Trimmer is a utility tool written in Go that trims adapter sequences from small RNA reads. The application is designed to a handle compressed fastq file (.fastq.gz) as an input and produces a compressed fastq file as an output.

## Features

- Adapter trimming for small RNA reads
- Optional base trimming from both 5' and 3' ends, useful for NEXTflex library preparation
- Input validation for adapter sequences
- Efficient processing using Go's concurrency capabilities

## Requirements

- Go version 1.x of installing from source
- Access to the terminal/command line

## Installation from source (Go required)

1. Clone the repository to your local machine using `git clone https://github.com/sfletc/scramTrimmer.git.
2. Navigate to the project directory with `cd scramTrimmer`.
3. Build the project using `go build`.

## Usage

```
./scramTrimmer -i inputfile.fastq.gz -o outputfile.fastq.gz -a adapter_sequence
```

**Parameters:**

- `-i`: Input file (required)
- `-o`: Output file (required)
- `-a`: Adapter sequence (required)
- `-minLen`: Minimum length of read after trimming (default 18)
- `-trim5`: 5' trim length (default 0)
- `-trim3`: 3' trim length after adapter removal (default 0)
- `-min5Match`: Minimum match length at 5' end (default 8)
- `-maxError`: Maximum mean error rate (default 0.1)

## Contribution

Contributions are welcome! Please make a pull request and we will review your code.

## License

BSD 3-Clause License

## Contact

- Stephen Fletcher: [s.fletcher@uq.edu.au](mailto:s.fletcher@uq.edu.au)
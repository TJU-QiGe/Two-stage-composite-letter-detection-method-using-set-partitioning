![cpp version](https://img.shields.io/badge/c%2B%2B-11-orange)
![shell version](https://img.shields.io/badge/shell-Bash-yellowgreen)
![license](https://img.shields.io/badge/license-MIT-8A2BE2)

[![Static Badge](https://img.shields.io/badge/github-TJU--QiGe-54af7d)](https://github.com/TJU-QiGe/Two-stage-composite-letter-detection-method-using-set-partitioning-for-DNA-storage-system)

<h1 id="TJU-QiGe">Two-stage composite letter detection method using set partitioning for DNA storage system</h1>

## Table of Contents

- [About the Project](#about-the-project)
- [Project Structure](#project-structure)
- [Encoding Pipeline](#encoding-pipeline)
- [Decoding Pipeline](#decoding-pipeline)
- [Files](#files)
- [Core Programs](#core-programs)
- [Example of Usage](#example-of-usage)
- [Dependencies](#dependencies)
- [License](#license)

## About the Project

DNA data storage, with its high-density potential and remarkable stability, is emerging as a strong candidate for future data storage. Composite DNA letters, a representation composed of a mixture of all four natural bases in a predetermined ratio, provide a approach to increase the logical density of DNA storage system. Here, we propose a two-stage composite letter detection method for the composite DNA storage system, achieving low-coverage data readout.

This project provides a complete DNA storage system implementation including both **encoding** and **decoding** pipelines:

- **Encoding Pipeline**: Converts user data files to composite DNA strands using Reed-Solomon error correction and composite letter mapping (for example, ATGCRYMK alphabet)
- **Decoding Pipeline**: Reconstructs original data from sequencing reads using two-stage composite letter detection and error correction

The system achieves a logical density of 2.5 bits per synthesis cycle using eight composite letters and demonstrates successful data readout from NGS sequencing data.

## Project Structure

```
src_upload_RS/
├── 126composite_strands_encoder/          # Encoding pipeline
│   ├── DNA_Encoding/                      # Encoding programs
│   │   ├── txt2bits.cpp                   # Text to binary conversion
│   │   ├── RS_1575_1890_encoder.cpp       # Reed-Solomon encoding with composite mapping
│   │   └── index_8nt_left.txt             # Left-end index sequences
│   │   └── index_8nt_right.txt            # Right-end index sequences
│   ├── Encoding_result/                   # Encoding output directory
│   └── User_Data/                         # Input data directory
│       └── Tang_Poems.txt                 # Example text file for encoding
├── Letter_detection_pipeline/             # Decoding pipeline
│   ├── full_readout_pipeline.sh           # Complete decoding pipeline script
│   ├── src/                               # Decoding programs source code
│   └── configureFiles/                    # Reference and configuration files
└── Demo_FASTQ_eight_letters/              # Example sequencing data
    └── Example-paired-end-assembled_100x.fastq
```

## Encoding Pipeline

The encoding pipeline converts a 2,355-byte Tang poem text file into composite DNA strands through a systematic process. The digital data are first converted into a 18,900-bit binary sequence, which is then encoded using a Reed–Solomon code RS(1890, 1575) over *GF*(2<sup>12</sup>) to introduce redundancy for error correction, producing a total of 22,680 bits. The encoded bitstream is divided into 7,560 chunks of 3 bits and mapped to octal symbols, which are subsequently assembled into 126 payload sequences of 60 letters each. These payloads correspond to 126 composite DNA strands composed of both natural bases and composite letters from the eight-letter alphabet {A, T, G, C, R, Y, M, K}.

### Encoding Steps

1. **Text to Binary Conversion**
   - Converts the input text file into a binary bit sequence
   - Pads the sequence to the required length (18,900 bits)

2. **Data Scrambling**
   - Applies XOR-based scrambling using a pseudo-random sequence
   - Improves data security and ensures uniform error distribution

3. **Reed-Solomon Encoding**
   - Performs RS(1890, 1575) over *GF*(2<sup>12</sup>)
   - Provides error correction with a coding rate of 0.833

4. **Composite Letter Mapping**
   - Maps each 3-bit pattern to one of eight composite letters {A, T, G, C, R, Y, M, K}
   - Generates 126 composite DNA strands, each 60 letters

5. **Index and Primer Addition**
   - Appends 8-nucleotide double-end indices for strand identification
   - Incorporates PCR primer sequences for amplification and sequencing

### Encoding Programs

| **Program** | **Function** | **Input** | **Output** |
|-------------|-------------|-----------|------------|
| `txt2bits` | Text to binary conversion | Text file, expected bits | Binary bit sequence |
| `RS_1575_1890_encoder` | Reed-Solomon encoding and composite mapping over *GF*(2<sup>12</sup>) | Binary data, scramble sequence, indices, output directory | Composite DNA strands with primers and indices |

### Composite Letter Alphabet

The system uses 8 composite letters with the following 3-bit mapping:

| **Binary** | **Composite Letter** | **Composition** |
|:----------:|:-------------------:|:---------------:|
| 000 | A | Pure Adenine |
| 001 | T | Pure Thymine |
| 010 | G | Pure Guanine |
| 011 | C | Pure Cytosine |
| 100 | R | A+G mixture  |
| 101 | Y | C+T mixture  |
| 110 | M | A+C mixture  |
| 111 | K | G+T mixture  |

## Decoding Pipeline

The decoding pipeline reconstructs original data from sequencing reads through five sequential steps:

### Decoding Steps

1. **Random Downsampling**
   - Randomly select reads from large sequencing datasets to achieve fixed coverage

2. **Primer Identification**
   - Identify and remove primer sequences through paired-end primer alignment

3. **Index Identification**
   - Cluster and group sequencing reads using double-end indices

4. **Two-stage Composite Letter Detection**
   - Apply set partitioning-based composite letter detection algorithm
   - Perform consensus letter detection through entropy calculation and maximum likelihood estimation

5. **Reed-Solomon Decoding**
   - Perform error correction and deletion recovery using Reed-Solomon error correction codes

## Files

### Table 1. Reference data and sequencing data used in our study.

| **Files**                             | **Storage Location**      | **Description**                                                     |
| ------------------------------------------- | ------------------------------- | ------------------------------------------------------------------------- |
| **Encoding Input Files** |  |  |
| `Tang_Poems.txt`                          | `./126composite_strands_encoder/User_Data/` | Example text file for encoding demonstration |
| `watermark_18900bit.txt`                  | `./configureFiles`            | Watermark sequence for data scrambling during encoding |
| `index_8nt_left.txt`                      | `./126composite_strands_encoder/DNA_Encoding/` | Left-end index sequences for strand identification |
| `index_8nt_right.txt`                     | `./126composite_strands_encoder/DNA_Encoding/` | Right-end index sequences for strand identification |
| **Decoding Reference Files** |  |  |
| `Tang_poem.txt`                           | `./configureFiles`            | The original text file stored in composite DNA strands |
| `composite_strand_ref.txt`                | `./configureFiles`            | Encoded composite payload sequences, length: 60 letters, Number: 126 |
| `index_8nt.txt`                           | `./configureFiles`            | The designed 5'-end index sequences for payload identification |
| `index_8nt_interleaved.txt`               | `./configureFiles`            | The designed 3'-end index sequences for payload identification |
| **Sequencing Data** |  |  |
| `Example-paired-end-assembled_100x.fastq` | `../Demo_FASTQ_eight_letters` | Paired-end assembled sequencing data based on NGS reads at 100x coverage |

## Core Programs

### Table 2. Encoding programs and their functions.

| **Program** | **Input** | **Output** | **Description** |
|-------------|-----------|------------|-----------------|
| `txt2bits` | 1. Input text file<br>2. Output file path<br>3. Expected total bits (18900) | Binary bit sequence file | Converts text files to binary representation, padding with zeros if necessary |
| `RS_1575_1890_encoder` | 1. Input data file (18900 bits)<br>2. Scramble sequence file<br>3. Output directory<br>4. Left index file<br>5. Right index file | 1. `126composite_strands_116letters.txt`<br>2. `original_bits_18900.txt`<br>3. `rs_encoded_codeword.txt`<br>4. `rscode_bits.txt`<br>5. `composite_payload_ref.txt` | Performs Reed-Solomon encoding, composite letter mapping, and generates complete DNA strands with primers and indices |

### Table 3. Decoding programs and their functions.

| **Program**         | **Input**                                                                                                                     | **Output**                                                                                                                                        | **Description**                                                                                                                             |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `seqtk`                 | 1. `Example-paired-end-assembled_100x.fastq` <br> 2. `Coverage=20`                                                        | `sub_sample.fastq`                                                                                                                                    | **Function 1**: Randomly downsamples sequencing reads from 100x to specified coverage using seqtk tool for reduced data processing volume.    |
| `primer_identification` | 1. `sub_sample.fastq` <br> 2. Forward primer: `TCACCATCCACTCTAAACAC` <br> 3. Reverse primer: `ATGAGTGGAGGTGTAAAGTG`    | `read_without_primer.fasta`                                                                                                                           | **Function 2**: Aligns forward and reverse primers with sequencing reads for primer identification and removal.                                 |
| `groupingByIndex`       | 1. `index_8nt.txt` <br> 2. `index_8nt_interleaved.txt` <br> 3. `read_without_primer.fasta` <br> 4. `Threshold=2` | 1. Grouped reads: `payload.txt` <br> 2. Match information: `matched_information_record.txt` <br> 3. Address frequency: `addrFrequency.txt` | **Function 3**: Performs double-end index identification to extract the payload region.              |
| `letter_detection_core` | 1. `payload.txt` <br> 2. `alphabet=1` (ATGCRYMK: Eight letters used in the composite DNA pool)                                                            | `consensus.txt`                                                                                                                                       | **Function 4**: Executes the two-stage composite letter detection algorithm using entropy-based partitioning and maximum likelihood estimation. |
| `data_analysis`         | 1. `consensus.txt` <br> 2. `payload.txt` <br> 3. `composite_strand_ref.txt` <br> 4. `basetype=ATGCRYMK`          | 1. `Different_letter_errors.txt` <br> 2. `Errors_per_strand.txt` <br> 3. `Total_errors.txt` <br> 4. `Entropy-*.txt`                  | **Function 5 (Optional)**: Analyzes detection accuracy by comparing consensus sequences with reference, calculating error rates and entropy distributions (used for performance evaluation). |
| `RS_decoder`            | 1. `consensus.txt` <br> 2. `alphabet=1`                                                                                      | `decoded_information_bits.txt`                                                                                                                        | **Function 6**: Performs Reed-Solomon error correction decoding on the consensus sequences to recover error-free information bits.              |
| `recovery_poem`         | 1. `decoded_information_bits.txt` <br> 2. `watermark_18900bit.txt`                                                           | `Tang_poem.txt`                                                                                                                                       | **Function 7**: Converts decoded binary data back to original poem file.             |

## Example of Usage

### Running the Encoding Pipeline

The encoding pipeline has been integrated into a single script for easy execution. Simply configure the parameters and run with one command:

```bash
cd 126composite_strands_encoder/DNA_Encoding/
bash run_DNA_Encoding.sh
```

The integrated script automatically handles all encoding steps including compilation, file preparation, and execution. Below are the detailed steps that the script performs:

#### Step 1: Prepare Input Data
Place your text file in the `User_Data` directory:
```bash
cd 126composite_strands_encoder/User_Data/
# Place your text file (e.g., Tang_Poems.txt) here
```

#### Step 2: Convert Text to Binary
```bash
cd ../DNA_Encoding/
g++ -o txt2bits txt2bits.cpp
./txt2bits ../User_Data/Tang_Poems.txt input_data_bits.txt 18900
```

#### Step 3: Perform Reed-Solomon Encoding and Composite Mapping
```bash
g++ -o RS_encoder RS_1575_1890_encoder.cpp -lm
mkdir -p ../Encoding_result
./RS_encoder input_data_bits.txt ../../configureFiles/watermark_18900bit.txt ../Encoding_result index_8nt_left.txt index_8nt_right.txt
```

#### Expected Encoding Output
```
126composite_strands_encoder/Encoding_result/
├── 126composite_strands_116letters.txt    # Complete DNA strands with primers and indices
├── original_bits_18900.txt                # Original binary data matrix
├── rs_encoded_codeword.txt                # Reed-Solomon codewords
├── rscode_bits.txt                        # RS-encoded binary data
└── composite_payload_ref.txt              # Composite letter payload sequences
```

### Running the Decoding Pipeline

Execute the integrated pipeline script:

```bash
cd Letter_detection_pipeline
bash full_readout_pipeline.sh
```

The pipeline will automatically execute all five steps and provide detailed progress information:

#### Step 1: Random Downsampling

- **Input**: `../Demo_FASTQ_eight_letters/Example-paired-end-assembled_100x.fastq`
- **Process**: Uses `seqtk` to downsample from 100x to 20x coverage (2,520 reads)
- **Output**: `../Eight_letter_result/20/downsampling/sub_sample.fastq`

#### Step 2: Primer Identification and Removal

- **Input**: Downsampled FASTQ file and primer sequences
- **Process**: `primer_identification` program aligns and removes primers
- **Output**: `../Eight_letter_result/20/downsampling/read_without_primer.fasta`

#### Step 3: Index-based Read Grouping

- **Input**: Primer-trimed reads and double-end index sequences
- **Process**: `groupingByIndex` identifies double-end indices and extracts payload regions
- **Output**: Grouped payload sequences in `../Eight_letter_result/20/group_reads/payload.txt`
- **Additional outputs**:
  - `matched_information_record.txt`: Statistical record file of index identification success rate
  - `addrFrequency.txt`: Read copy number distribution per composite strand

#### Step 4: Composite Letter Detection and Analysis

- **Input**: Payload reads grouped by double-end index
- **Process**:
  - `letter_detection_core` performs two-stage detection using set partitioning based on entropy of observed base frequency
  - `data_analysis` computes per-letter and per-strand accuracy and error statistics
- **Output**:
  - `consensus.txt`: consensus sequences per strand
  - Error analysis files: `Different_letter_errors.txt`, `Errors_per_strand.txt`, `Total_errors.txt`
  - Entropy distribution files: `Entropy-*.txt` (one file per base)

#### Step 5: Reed-Solomon Decoding and Data Recovery

- **Input**: Consensus sequences with errors
- **Process**:
  - `RS_decoder` corrects residual errors via Reed–Solomon block decoding
  - `recovery_poem` converts the corrected bitstream back to the original text
- **Output**:
  - `decoded_information_bits.txt`: Error-corrected information bits
  - `Tang_poem.txt`: reconstructed source text

### Complete System Workflow

For a complete encode-decode cycle:

```bash
# 1. Encoding: Text → Composite DNA Strands
cd 126composite_strands_encoder/DNA_Encoding/
g++ -o txt2bits txt2bits.cpp
g++ -o RS_encoder RS_1575_1890_encoder.cpp -lm
mkdir -p ../Encoding_result

./txt2bits ../User_Data/Tang_Poems.txt input_data_bits.txt 18900
./RS_encoder input_data_bits.txt ../../configureFiles/watermark_18900bit.txt ../Encoding_result index_8nt_left.txt index_8nt_right.txt

# 2. Simulation: DNA Strands → Sequencing Data (external NGS simulation)

# 3. Decoding: Sequencing Data → Original Text
cd ../../Letter_detection_pipeline/
bash full_readout_pipeline.sh
```

### Configuration Options

Both pipelines can be customized:

**Encoding Configuration:**
- Input text file size and content
- Reed-Solomon parameters (currently RS(1890,1575))
- Composite letter alphabet (currently ATGCRYMK)
- Index sequences for strand identification
- Primer sequences for PCR amplification

**Decoding Configuration:**
- `COVERAGE`: Target sequencing coverage (default: 20)
- `THRESHOLD`: Index alignment threshold (default: 2)
- `BASE_TYPES`: Composite letter alphabet (default: "ATGCRYMK")
- `ALPHABET`: Composite base complexity (default: 1)

### Performance Monitoring

The pipeline provides comprehensive execution monitoring:

- Real-time progress updates for each step
- Compilation status for all programs
- File generation confirmations
- Execution time measurements (CPU and wall time)
- Memory usage statistics
- Final summary with complete file listing

## Dependencies

The following external tools are required for full functionality:

- **seqtk**: For subsampling and processing FASTQ files
  - URL: [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)
  
- **edlib**: For sequence alignment operations
  - URL: [https://github.com/Martinsos/edlib](https://github.com/Martinsos/edlib)

**System Requirements:**
- Linux operating system (tested on Ubuntu 18.04+)
- C/C++ compiler supporting C++11 standard
- Bash shell environment
- Sufficient memory for processing large sequencing files

**Installation:**
```bash
# Install required tools
git clone https://github.com/lh3/seqtk.git
cd seqtk && make

git clone https://github.com/Martinsos/edlib.git
cd edlib && make
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

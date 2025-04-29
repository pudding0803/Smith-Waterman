# Smith-Waterman

## Initialize Submodules

```bash
git submodule update --init --recursive
```

## Build and Run

```bash
mkdir build
cd build
cmake ..
cmake --build .
./main
```

## Example Data

> Example inputs are located in the `/data` directory.

- From [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
  - `ACTB.fa`
  - `BCL2.fa`
  - `GAPDH.fa`
  - `TP53.fa`

## Supported Modes

- `SmithWaterman::Mode`
  - `Naive`
  - `SIMD`

## Usage

```cpp
// main.cpp
auto report = SmithWaterman::run(fasta1, fasta2, mode);
```

- Structure of `AlignmentReport`
  - `alignmentResult`
  - `durationMs`
  - `maxScore`


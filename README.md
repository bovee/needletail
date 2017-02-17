# OLDER VERSION

This is an older, unmaintained version of Needletail. Please find the current version [here](https://github.com/onecodex/needletail).

# Needletail

Needletail is a MIT licensed, minimal-copying FASTA/FASTQ parser and k-mer processing library.

The goal is to write a fast *and* well-tested set of functions that more-specialized bioinformatics programs can use. Needletail's goal is to be as fast as the `readfq` C library at parsing FASTX files and much (i.e. 25 times) faster than equivalent Python implementations at k-mer counting.

# Example

```rust
extern crate needletail;
ues std::env;
use needletail::{fastx};

fn main() {
  let filename: String = env::args().nth(1).unwrap();

  let mut n_bases = 0;
  fastx::fastx_file(&filename[..], |seq| {
    // seq.0 is the name of the record
    // seq.1 is the base sequence
    n_bases += seq.1.len();
    // seq.2 is an optional quality score
  });
  // the below number will include line endings, so it's only really
  // valid for FASTQs
  println!("There are {} bases in your file.", n_bases);
}
```

# Installation

Needletail requires `rust` and `cargo` to be installed.

```shell
git clone https://github.com/bovee/needletail
cargo test  # to run tests
```

# Getting Help

Questions are best directed as GitHub issues.

Hopefully I'll compile the documentation and put it up as a webpage soon too.

# Contributing

Please do! I'm happy to discuss/mentor possible additions and/or accept pull requests.

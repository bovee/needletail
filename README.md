# Needletail

Needletail is a MIT licensed, minimal-copying FASTA/FASTQ parser and k-mer processing library.

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
  println!("{}", n_bases);
}
```


# Installation

Needletail requires `rust` and `cargo` to be installed.

```shell
git clone https://github.com/bovee/needletail
cargo test  # to run tests
```

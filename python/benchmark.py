#!/usr/bin/env python
"""
Unscientific benchmarking of this versus the --release rust
implementation below using the %timeit Ipython magic (times in sec)

    n_kmers,  py_runtime, rust_runtime
    736870,   1.34,       0.0965
    3988647,  6.43,       0.417
    99224925, 161,        10.2

Both give identical counts on the files tested (and printing kmers out
and diff'ing the results gives no difference)

```rust
extern crate needletail;
use std::string::String;
use needletail::{fastx, kmer};

fn main() {
    let mut n_total = 0;
    let mut n_canonical = 0;
    fastx::fastx_file("../../../Downloads/20160428.curated.loci.fasta", |seq| {
        for k in kmer::Kmer::new(seq.1, 4, true) {
            let l = kmer::canonical(k).into_owned();
            // println!("{:}", String::from_utf8_lossy(&l));
            if l == &b"CAGC"[..] {
                n_canonical += 1;
            }
            n_total += 1;
        }
    });
    println!("{:?} {:?}", n_total, n_canonical);
}
```

"""
from __future__ import print_function
from Bio.SeqIO import parse
from Bio.Seq import reverse_complement


def slid_win(seq, size=4, overlapping=True):
    """Returns a sliding window along self.seq."""
    itr = iter(seq)
    if overlapping:
        buf = ''
        for _ in range(size):
            buf += next(itr)
        for l in itr:
            yield buf
            buf = buf[1:] + l
        yield buf
    else:
        buf = ''
        for l in itr:
            buf += l
            if len(buf) == size:
                yield buf
                buf = ''


n_total = 0
n_canonical = 0

for s in parse('test_file.fa', 'fasta'):
    uppercase_seq = str(s.upper().seq)
    for kmer in slid_win(uppercase_seq, 4):
        canonical = min(kmer, reverse_complement(kmer))
        if canonical == 'CAGC':
            n_canonical += 1
        n_total += 1
        # print(canonical)

print(n_total, n_canonical)

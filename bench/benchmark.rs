extern crate needletail;

use std::env;
use needletail::{fastx, kmer};


fn count_bases(filename: String) {
    let mut n_bases = 0;
    fastx::fastx_file(&filename[..], |seq| {
        n_bases += seq.1.len();
    });
    println!("{:?}", n_bases);
}


fn count_canonical_kmers(filename: String) {
    let mut n_total = 0;
    let mut n_canonical = 0;
    fastx::fastx_file(&filename[..], |seq| {
        for k in kmer::Kmer::new(seq.1, 31, true) {
            let l = kmer::canonical(k).into_owned();
            if l == &b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"[..] {
                n_canonical += 1;
            }
            n_total += 1;
        }
    });
    println!("{:?} {:?}", n_total, n_canonical);
}

fn main() {
    let filename: String = env::args().nth(1).unwrap();

    count_bases(filename);
}

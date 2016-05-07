extern crate needletail;
use std::string::String;
use needletail::{fastx, kmer};


fn main() {
    let mut n_total = 0;
    let mut n_canonical = 0;
    fastx::fastx_file("../../../Downloads/20160428.curated.loci.fasta", |seq| {
        for k in kmer::Kmer::new(seq.1, 4, true) {
            // let l = kmer::normalize(kmer::canonical(k)).into_owned();
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

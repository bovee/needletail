extern crate needletail;

use std::borrow::Cow;
use std::env;
use needletail::{fastx, kmer};

// fn main() {
//     let filename: String = env::args().nth(1).unwrap();
//
//     let mut n_total = 0;
//     let mut n_canonical = 0;
//     fastx::fastx_file(&filename[..], |seq| {
//         // let normalized_seq = kmer::normalize(Cow::Borrowed(&seq.1), true);
//         for k in seq.1.windows(31) {
//             let l = kmer::canonical(k).into_owned();
//             //            if l == &b"CAGC"[..] {
//             //                n_canonical += 1;
//             //            }
//             n_total += 1;
//         }
//     });
//     println!("{:?} {:?}", n_total, n_canonical);
// }

fn main() {
    let filename: String = env::args().nth(1).unwrap();

    let mut n_bases = 0;
    fastx::fastx_file(&filename[..], |seq| {
        n_bases += seq.1.len();
    });
    println!("{:?}", n_bases);
}

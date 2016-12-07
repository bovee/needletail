#[macro_use]
extern crate bencher;
extern crate needletail;

use needletail::{fastx, kmer};
use bencher::Bencher;

fn check_kmer_speed(bench: &mut Bencher) {
    let filename = String::from("tests/data/28S.fasta");

    bench.iter(|| {
        let mut n_total = 0;
        let mut n_canonical = 0;
        fastx::fastx_file(&filename[..], |seq| {
            // let normalized_seq = kmer::normalize(Cow::Borrowed(&seq.1), true);
            for k in seq.1.windows(31) {
                let l = kmer::canonical(k).into_owned();
                if l == k {
                    n_canonical += 1;
                }
                n_total += 1;
            }
        });
        assert_eq!(213703, n_total);
        assert_eq!(108521, n_canonical);
    })
}

benchmark_group!(benches, check_kmer_speed);
benchmark_main!(benches);

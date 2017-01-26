use bitkmer::BitNuclKmer;
use kmer::{is_good_base, complement};


pub struct Seq<'a> {
    seq: &'a [u8],
    rev_seq: Vec<u8>,
}

impl<'a> Seq<'a> {
    pub fn new(seq: &'a [u8]) -> Self {
        let rev_seq = seq.iter().rev().map(|n| complement(n)).collect();
        Seq {
            seq: seq,
            rev_seq: rev_seq,
        }
    }

    // pub fn kmers(&'a self, k: u8) -> Windows<&[u8]> {
    //     self.seq.windows(k as usize)
    // }

    pub fn valid_kmers(&'a self, k: u8) -> NuclKmer<'a> {
        NuclKmer::new(self.seq, self.seq, k)
    }

    pub fn canonical_kmers(&'a self, k: u8) -> NuclKmer<'a> {
        NuclKmer::new(self.seq, &self.rev_seq, k)
    }

    pub fn bit_kmers(&'a self, k: u8) -> BitNuclKmer<'a> {
        BitNuclKmer::new(self.seq, k, false)
    }

    pub fn canonical_bit_kmers(&'a self, k: u8) -> BitNuclKmer<'a> {
        BitNuclKmer::new(self.seq, k, true)
    }
}

pub struct NuclKmer<'a> {
    k: u8,
    start_pos: usize,
    buffer: &'a [u8],
    rc_buffer: &'a [u8],
}

fn update_position(start_pos: &mut usize, k: u8, buffer: &[u8], initial: bool) -> bool {
    // check if we have enough "physical" space for one more kmer
    if *start_pos + k as usize > buffer.len() {
        return false;
    }

    let mut kmer_len = (k - 1) as usize;
    let mut stop_len = k as usize;
    if initial {
        kmer_len = 0;
        stop_len = (k - 1) as usize;
    }

    while kmer_len < stop_len {
        if is_good_base(buffer[*start_pos + kmer_len]) {
            kmer_len += 1;
        } else {
            kmer_len = 0;
            *start_pos += kmer_len + 1;
            if *start_pos + k as usize > buffer.len() {
                return false;
            }
        }
    }
    true
}

impl<'a> NuclKmer<'a> {
    //! A kmer-izer for a nucleotide/amino acid sequence; returning slices to the original data
    pub fn new(buffer: &'a [u8], rc_buffer: &'a [u8], k: u8) -> NuclKmer<'a> {
        let mut start_pos = 0;
        update_position(&mut start_pos, k, buffer, true);
        NuclKmer {
            k: k,
            start_pos: start_pos,
            buffer: buffer,
            rc_buffer: rc_buffer,
        }
    }
}

impl<'a> Iterator for NuclKmer<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<&'a [u8]> {
        if !update_position(&mut self.start_pos, self.k, self.buffer, false) {
            return None;
        }

        let result = &self.buffer[self.start_pos..self.start_pos + self.k as usize];
        let rc_result = &self.rc_buffer[self.start_pos..self.start_pos + self.k as usize];
        self.start_pos += 1;
        if result < rc_result {
            Some(result)
        } else {
            Some(rc_result)
        }
    }
}


#[test]
fn can_kmerize() {
    // test general function
    let mut i = 0;
    for k in Seq::new(b"AGCT").valid_kmers(1) {
        match i {
            0 => assert_eq!(k, &b"A"[..]),
            1 => assert_eq!(k, &b"G"[..]),
            2 => assert_eq!(k, &b"C"[..]),
            3 => assert_eq!(k, &b"T"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that we skip over N's
    i = 0;
    for k in Seq::new(b"ACNGT").valid_kmers(2) {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            1 => assert_eq!(k, &b"GT"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that we skip over N's and handle short kmers
    i = 0;
    for k in Seq::new(b"ACNG").valid_kmers(2) {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that the minimum length works
    for k in Seq::new(b"AC").valid_kmers(2) {
        assert_eq!(k, &b"AC"[..]);
    }
}

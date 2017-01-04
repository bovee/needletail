//! This module contains functions for kmerizing a longer sequence and various
//! utilities for dealing with these kmers.
//!
//!
use std::borrow::Cow;

#[inline]
fn complement(n: &u8) -> u8 {
    //! Returns the complementary base for a given IUPAC base code.
    //!
    //! Does not work for RNA sequences (maybe we should raise an error or something?)
    match *n as char {
        'a' => 't' as u8,
        'c' => 'g' as u8,
        'g' => 'c' as u8,
        't' => 'a' as u8,
        'A' => 'T' as u8,
        'C' => 'G' as u8,
        'G' => 'C' as u8,
        'T' => 'A' as u8,

        // IUPAC codes
        'r' => 'y' as u8,
        'y' => 'r' as u8,
        's' => 's' as u8,
        'w' => 'w' as u8,
        'k' => 'm' as u8,
        'm' => 'k' as u8,
        'b' => 'v' as u8,
        'd' => 'h' as u8,
        'h' => 'd' as u8,
        'v' => 'b' as u8,
        'R' => 'Y' as u8,
        'Y' => 'R' as u8,
        'S' => 'S' as u8,
        'W' => 'W' as u8,
        'K' => 'M' as u8,
        'M' => 'K' as u8,
        'B' => 'V' as u8,
        'D' => 'H' as u8,
        'H' => 'D' as u8,
        'V' => 'B' as u8,

        // anything else just pass through
        // 'u' | 'U' => panic!("Does not support complements of U"),
        x => x as u8,
    }
}

pub fn normalize<'a>(seq: &'a [u8], iupac: bool) -> Cow<'a, [u8]> {
    //! Transform a FASTX sequence into it's "normalized" form.
    //!
    //! The normalized form is:
    //!  - only AGCTN and possibly . (for gaps)
    //!  - lowercase versions of these are uppercased
    //!  - U is converted to T (make everything a DNA sequence)
    //!  - some other punctuation is converted to gaps
    //!  - IUPAC bases may be converted to N's depending on the parameter passed in
    //!  - everything else is considered a N
    let mut buf: Vec<u8> = Vec::with_capacity(seq.len());
    let mut original_was_bad = false;

    for n in seq.iter() {
        let mut good_char = false;
        buf.push(match (*n as char, iupac) {
            c @ ('A', _) | c @ ('C', _) | c @ ('G', _) | c @ ('T', _) | c @ ('N', _) |
            c @ ('.', _) => {
                good_char = true;
                c.0 as u8
            },
            ('a', _) => 'A' as u8,
            ('c', _) => 'C' as u8,
            ('g', _) => 'G' as u8,
            // normalize uridine to thymine
            ('t', _) | ('u', _) | ('U', _) => 'T' as u8,
            ('-', _) | ('~', _) | (' ', _) => '.' as u8,
            // logic for IUPAC bases (a little messy)
            c @ ('B', true) | c @ ('D', true) | c @ ('H', true) | c @ ('V', true) |
            c @ ('R', true) | c @ ('Y', true) | c @ ('S', true) | c @ ('W', true) |
            c @ ('K', true) | c @ ('M', true) => {
                good_char = true;
                c.0 as u8
            },
            ('b', true) => 'B' as u8,
            ('d', true) => 'D' as u8,
            ('h', true) => 'H' as u8,
            ('v', true) => 'V' as u8,
            ('r', true) => 'R' as u8,
            ('y', true) => 'Y' as u8,
            ('s', true) => 'S' as u8,
            ('w', true) => 'W' as u8,
            ('k', true) => 'K' as u8,
            ('m', true) => 'M' as u8,
            _ => 'N' as u8,
        });
        if !good_char {
            original_was_bad = true;
        }
    }
    match original_was_bad {
        true => Cow::Owned(buf),
        false => Cow::Borrowed(seq),
    }
}

pub fn canonical<'a>(seq: &'a [u8]) -> Cow<'a, [u8]> {
    //! Taking in a sequence string, return the canonical form of the sequence
    //! (e.g. the lexigraphically lowest of either the original sequence or its
    //! reverse complement)
    let mut buf: Vec<u8> = Vec::with_capacity(seq.len());
    // enough just keeps our comparisons from happening after they need to
    let mut enough = false;
    let mut original_was_canonical = false;

    // loop through the kmer and its reverse complement simultaneously
    for (rn, n) in seq.iter().rev().map(|n| complement(n)).zip(seq.iter()) {
        buf.push(rn);
        if !enough && n < &rn {
            original_was_canonical = true;
            break;
        } else if !enough && &rn < n {
            enough = true;
        }
        // unstated if branch: if rn == n, keep comparing
    }
    match (original_was_canonical, enough) {
        (true, true) => panic!("Bug: should never set original_was_canonical if enough == true"),
        (true, false) => Cow::Borrowed(seq),
        (false, true) => Cow::Owned(buf),
        // the sequences were completely equal, return the ref
        (false, false) => Cow::Borrowed(seq),
    }
}

#[test]
fn can_canonicalize() {
    assert!(canonical(b"A") == Cow::Borrowed(b"A"));
    assert!(canonical(b"T") == Cow::Owned::<[u8]>(b"A".to_vec()));
    assert!(canonical(b"AAGT") == Cow::Borrowed(b"AAGT"));
    assert!(canonical(b"ACTT") == Cow::Owned::<[u8]>(b"AAGT".to_vec()));
    assert!(canonical(b"GC") == Cow::Borrowed(b"GC"));
}

/// Find the lexigraphically smallest substring of `seq` of length `length`
///
/// There's probably a faster algorithm for this somewhere...
pub fn minimizer<'a>(seq: &'a [u8], length: usize) -> Cow<'a, [u8]> {
    let reverse_complement: Vec<u8> = seq.iter().rev().map(|n| complement(n)).collect();
    let mut minmer = Cow::Borrowed(&seq[..length]);

    for (kmer, rc_kmer) in seq.windows(length).zip(reverse_complement.windows(length)) {
        if kmer < &minmer[..] {
            minmer = Cow::Borrowed(kmer);
        }
        if rc_kmer < &minmer[..] {
            minmer = Cow::Owned(rc_kmer.to_vec());
        }
    }
    minmer
}

#[test]
fn can_minimize() {
    let minmer = minimizer(&b"ATTTCG"[..], 3);
    assert_eq!(&minmer[..], b"AAA");
}

// TODO
// pub fn skip_n<'a, T>(iter: T) -> T where T: Iterator<Item=&'a [u8]> {
//    iter.filter(|kmer| kmer.contains(&('N' as u8)) || kmer.contains(&('n' as u8)))
// }

fn is_good_base(chr: u8) -> bool {
    match chr as char {
        'a' | 'c' | 'g' | 't' | 'A' | 'C' | 'G' | 'T' => true,
        _ => false,
    }
}

fn has_no_n<'a>(seq: &'a [u8]) -> bool {
    //! Determines if a sequence has any non-primary four bases
    //! characters in it
    seq.iter().all(|n| is_good_base(*n))
}

#[test]
fn can_detect_no_n() {
    assert!(has_no_n(b"AAGT"));
    assert!(!has_no_n(b"NAGT"));
}

pub struct NuclKmer<'a> {
    k: u8,
    start_pos: usize,
    buffer: &'a [u8],
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
    pub fn new(buffer: &'a [u8], k: u8) -> NuclKmer<'a> {
        let mut start_pos = 0;
        update_position(&mut start_pos, k, buffer, true);
        NuclKmer {
            k: k,
            start_pos: start_pos,
            buffer: buffer,
        }
    }
}

impl<'a> Iterator for NuclKmer<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<&'a [u8]> {
        if !update_position(&mut self.start_pos, self.k, self.buffer, false) {
            return None;
        }

        let result = Some(&self.buffer[self.start_pos..self.start_pos + self.k as usize]);
        self.start_pos += 1;
        result
    }
}


#[test]
fn can_kmerize() {
    // test general function
    let mut i = 0;
    for k in NuclKmer::new(b"AGCT", 1) {
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
    for k in NuclKmer::new(b"ACNGT", 2) {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            1 => assert_eq!(k, &b"GT"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that we skip over N's and handle short kmers
    i = 0;
    for k in NuclKmer::new(b"ACNG", 2) {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that the minimum length works
    for k in NuclKmer::new(b"AC", 2) {
        assert_eq!(k, &b"AC"[..]);
    }
}

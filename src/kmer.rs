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

pub fn normalize<'a>(seq: Cow<'a, [u8]>, iupac: bool) -> Cow<'a, [u8]> {
    //! Transform a FASTX sequence into it's "normalized" form.
    //!
    //! The normalized form is:
    //!  - only AGCTN and possibly . (for gaps)
    //!  - lowercase versions of these are uppercased
    //!  - U is converted to T (make everything a DNA sequence)
    //!  - some other punctuation is converted to gaps
    //!  - whitespace is removed (except spaces, which are considered gaps)
    //!  - IUPAC bases may be converted to N's depending on the parameter passed in
    //!  - everything else is considered a N
    let mut buf: Vec<u8> = Vec::with_capacity(seq.len());
    let mut original_was_bad = false;

    for n in seq.iter() {
        let mut good_char = false;
        buf.push(match (*n as char, iupac) {
            c @ ('A', _) | c @ ('C', _) | c @ ('G', _) | c @ ('T', _) | c @ ('N', _) | c @ ('.', _) => {
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
            ('\r', _) | ('\n', _) | ('\t', _) => {
                original_was_bad = true;
                // continue to prevent pushing whitespace into the normalized form
                continue;
            },
            _ => 'N' as u8,
        });
        if !good_char {
            original_was_bad = true;
        }
    }
    match original_was_bad {
        true => Cow::Owned(buf),
        false => seq,
    }
}

pub fn canonical<'a>(seq: Cow<'a, [u8]>) -> Cow<'a, [u8]> {
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
        (true, false) => seq,
        (false, true) => Cow::Owned(buf),
        (false, false) => seq  // the sequences were completely equal, return the ref
    }
}

#[test]
fn can_canonicalize() {
    // TODO: figure out a way to compare a Cow::Owned?
    //assert!(canonical(Cow::Borrowed(b"A")) == Cow::Owned::<[u8]>(b"T".to_vec()));
    assert!(canonical(Cow::Borrowed(b"AAGT")) == Cow::Borrowed(b"AAGT"));
    assert!(canonical(Cow::Borrowed(b"GC")) == Cow::Borrowed(b"GC"));
}

// TODO
//pub fn skip_n<'a, T>(iter: T) -> T where T: Iterator<Item=&'a [u8]> {
//    iter.filter(|kmer| kmer.contains(&('N' as u8)) || kmer.contains(&('n' as u8)))
//}

fn has_no_n(seq: Cow<[u8]>) -> bool {
    //! Determines if a sequence has any non-primary four bases
    //! characters in it
    seq.iter().all(|n| match *n as char {
        'a' | 'c' | 'g' | 't' | 'A' | 'C' | 'G' | 'T' => true,
        _ => false
    })
}

#[test]
fn can_detect_no_n() {
    assert!(has_no_n(Cow::Borrowed(b"AAGT")));
    assert!(!has_no_n(Cow::Borrowed(b"NAGT")));
}

pub struct Kmer<'a> {
    k: usize,
    multiline: bool,
    start_pos: Option<usize>,
    buffer: &'a [u8],
}

impl<'a> Kmer<'a> {
    //! A kmer-izer for a nucleotide/amino acid sequence; either returning slices to the
    //! original data or transparently (using CoW) removing whitespace
    pub fn new(slice: &'a [u8], k: usize, multiline: bool) -> Kmer<'a> {
        Kmer {
            k: k,
            multiline: multiline,
            start_pos: Some(0),
            buffer: slice,
        }
    }
}

impl<'a> Iterator for Kmer<'a> {
    type Item = Cow<'a, [u8]>;

    fn next(&mut self) -> Option<Cow<'a, [u8]>> {
        match self.start_pos {
            // the sequence is exhausted, return None as the Iterator sentinel value
            None => None,
            Some(pos) => {
                // check if we have enough "physical" space for one more kmer
                if pos > self.buffer.len() - self.k {
                    self.start_pos = None;
                    return None;
                }

                // if we're not in multiline mode, we can return whatever
                // (including whitespace)
                if !self.multiline {
                    self.start_pos = Some(pos + 1);
                    return Some(Cow::Borrowed(&self.buffer[pos..pos + self.k]));
                }
            
                let mut buf: Vec<u8> = Vec::with_capacity(self.k);
                let mut had_whitespace = false;
                let mut start_pos = pos;

                for n in self.buffer[pos..].iter() {
                    if n != &('\n' as u8) && n != &('\r' as u8) && n != &('\t' as u8) {
                        buf.push(*n);
                    } else if buf.len() <= 1 {
                        had_whitespace = true;
                        // bump up the start
                        start_pos += 1;
                    } else {
                        had_whitespace = true;
                    }
                    if buf.len() == self.k {
                        break;
                    }
                }
                
                // FIXME: this may not work for a "\r\n" because it's only going ahead one?
                // if the next position is going to start with whitespace, skip ahead
                if buf.len() != self.k {
                    self.start_pos = None;
                    None
                } else if !had_whitespace {
                    self.start_pos = Some(start_pos + 1);
                    Some(Cow::Borrowed(&self.buffer[pos..pos + self.k]))
                } else {
                    self.start_pos = Some(start_pos + 1);
                    Some(Cow::Owned(buf))
                }
            }
        }
    }
}


#[test]
fn can_kmerize() {
    // test general function
    let mut i = 0;
    for k in Kmer::new(b"AG\nCT", 1, true) {
        match i {
            0 => assert_eq!(k, &b"A"[..]),
            1 => assert_eq!(k, &b"G"[..]),
            2 => assert_eq!(k, &b"C"[..]),
            3 => assert_eq!(k, &b"T"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that we can join over newlines
    i = 0;
    for k in Kmer::new(b"A\nC\nG\nT\n", 2, true) {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            1 => assert_eq!(k, &b"CG"[..]),
            2 => assert_eq!(k, &b"GT"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that multiple non-space whitespaces work
    i = 0;
    for k in Kmer::new(b"A\r\nC \nG", 2, true) {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            1 => assert_eq!(k, &b"C "[..]),
            2 => assert_eq!(k, &b" G"[..]),
            _ => assert!(false),
        }
        i += 1;
    }
    
    // test that the minimum length works
    i = 0;
    for k in Kmer::new(b"A\r\nC", 2, true) {
        assert_eq!(k, &b"AC"[..]);
        assert!(i != 1);
        i += 1;
    }
}

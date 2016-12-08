pub type BitKmer = u64;

pub struct BitKmerIter<'a> {
    k: u8,
    pos: usize,
    cur_kmer: BitKmer,
    buffer: &'a [u8],
}

pub fn fast_extend_kmer(kmer: &BitKmer, k: &u8, new_char: &u8) -> BitKmer {
    // new_char must be in ACGTE (E gets counted as T unfortunately)
    assert!(new_char & 223 == 84 || new_char & 217 == 65);
    let new_kmer = (kmer << 2) + (((new_char >> 2) ^ (new_char >> 1)) & 3u8) as BitKmer;

    // mask out any overflowed bits
    new_kmer & (BitKmer::pow(2, (2 * k) as u32) - 1) as BitKmer
}

pub fn extend_kmer(kmer: &BitKmer, k: &u8, new_char: &u8) -> BitKmer {
    let new_kmer = (kmer << 2) + match new_char {
        &b'A' | &b'a' => 0 as BitKmer,
        &b'C' | &b'c' => 1 as BitKmer,
        &b'G' | &b'g' => 2 as BitKmer,
        &b'T' | &b't' => 3 as BitKmer,
        _ => panic!("Received a non ACGT nucleotide"),
    };

    // mask out any overflowed bits
    new_kmer & (BitKmer::pow(2, (2 * k) as u32) - 1) as BitKmer
}

impl<'a> BitKmerIter<'a> {
    pub fn new(slice: &'a [u8], k: u8) -> BitKmerIter<'a> {
        // TODO: raise an error if the slice is shorter than k
        let mut kmer = 0u64;
        if slice.len() >= k as usize {
            for i in 0..k - 1 {
                kmer = fast_extend_kmer(&kmer, &k, &slice[i as usize]);
            }
        }

        BitKmerIter {
            k: k,
            pos: (k - 2) as usize,
            cur_kmer: kmer,
            buffer: slice,
        }
    }
}

impl<'a> Iterator for BitKmerIter<'a> {
    type Item = BitKmer;

    fn next(&mut self) -> Option<BitKmer> {
        self.pos += 1;
        if self.pos > self.buffer.len() {
            return None;
        }
        self.cur_kmer = fast_extend_kmer(&self.cur_kmer, &self.k, &self.buffer[self.pos]);

        Some(self.cur_kmer)
    }
}

#[test]
fn test_iterator() {
    let seq = "ACGTA".as_bytes();
    let mut kmer_iter = BitKmerIter::new(seq, 3);
    assert_eq!(kmer_iter.next(), Some(6));
    assert_eq!(kmer_iter.next(), Some(27));
    assert_eq!(kmer_iter.next(), Some(44));
}

fn reverse_complement(kmer: BitKmer, k: u8) -> BitKmer {
    // FIXME: this is not going to work with BitKmers of u128 or u32
    // inspired from https://www.biostars.org/p/113640/
    let mut new_kmer = kmer;
    // reverse it
    new_kmer = (new_kmer >> 2 & 0x3333333333333333) | (new_kmer & 0x3333333333333333) << 2;
    new_kmer = (new_kmer >> 4 & 0x0F0F0F0F0F0F0F0F) | (new_kmer & 0x0F0F0F0F0F0F0F0F) << 4;
    new_kmer = (new_kmer >> 8 & 0x00FF00FF00FF00FF) | (new_kmer & 0x00FF00FF00FF00FF) << 8;
    new_kmer = (new_kmer >> 16 & 0x0000FFFF0000FFFF) | (new_kmer & 0x0000FFFF0000FFFF) << 16;
    new_kmer = (new_kmer >> 32 & 0x00000000FFFFFFFF) | (new_kmer & 0x00000000FFFFFFFF) << 32;
    // complement it
    new_kmer ^= 0xFFFFFFFFFFFFFFFF;
    // shift it to the right size
    new_kmer = new_kmer >> (2 * (32 - k));
    new_kmer
}

#[test]
fn test_reverse_complement() {
  assert_eq!(reverse_complement(0b000000, 3), 0b111111);
  assert_eq!(reverse_complement(0b111111, 3), 0b000000);
  assert_eq!(reverse_complement(0b00000000, 4), 0b11111111);
  assert_eq!(reverse_complement(0b00011011, 4), 0b00011011);
}

fn canonical(kmer: BitKmer, k: u8) -> BitKmer {
    let rc = reverse_complement(kmer, k);
    if kmer > rc {
        rc
    } else {
        kmer
    }
}

fn print_ikmer(kmer: BitKmer, k: u8) -> String {
    let mut new_kmer = kmer;
    let mut new_kmer_str = String::new();
    let offset = (k - 1) * 2;
    let bitmask = BitKmer::pow(2, (2 * k - 1) as u32) + BitKmer::pow(2, (2 * k - 2) as u32);

    for i in 0..k {
        let new_char = (new_kmer & bitmask) >> offset;
        new_kmer <<= 2;
        new_kmer_str.push(match new_char {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => panic!("Mathematical impossibility"),
        });
    }
    new_kmer_str
}

#[test]
fn test_print_ikmer() {
    assert_eq!(print_ikmer(1 as BitKmer, 1), String::from("C"));
    assert_eq!(print_ikmer(60 as BitKmer, 3), String::from("TTA"));
    assert_eq!(print_ikmer(0 as BitKmer, 3), String::from("AAA"));
}

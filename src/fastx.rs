//! This module contains functions for reading FASTA data from either within
//! memory or from files.
//!
//! # Design Note
//!
//! These functions are designed to take callbacks to process the FASTX records
//! they read. It would be nice to present a FASTX Iterator that downstream users
//! can use at some point, but this is only possible for in-memory data (using
//! MemProducer from nom) because otherwise the lifetime of each record is linked
//! to what part of the file we're reading and Rust doesn't support "streaming
//! iterators" otherwise. Maybe some day.
//!
//! See: https://github.com/emk/rust-streaming

use memchr::memchr;
use std::str;
use nom::{ErrorKind, IResult, Input, Move, Needed};
use nom::{Consumer, ConsumerState, FileProducer, MemProducer, Producer};

/// A generic FASTX record containing:
///   seq.0 - an id
///   seq.1 - the sequence itself
///   seq.2 - an optional set of quality scores for the sequence
pub type SeqRecord<'a> = (&'a str, Vec<u8>, Option<&'a [u8]>);

fn strip_whitespace(seq: &[u8]) -> Vec<u8> {
    seq.iter()
       .filter(|&n| *n != b'\r' && *n != b'\n')
       .map(|n| n.clone())
       .collect::<Vec<u8>>()
}

fn memchr2(b1: u8, b2: u8, seq: &[u8]) -> Option<usize> {
    // unlike memchr::memchr2, this looks for the bytes in sequence (not either/or)
    let mut pos = 0;
    while true {
        match memchr(b1, &seq[pos..]) {
            None => return None,
            Some(match_pos) => {
                if pos + match_pos + 1 == seq.len() {
                    return None;
                } else if seq[pos + match_pos + 1] == b2 {
                    return Some(pos + match_pos);
                } else {
                    pos += match_pos + 1;
                }
            },
        }
    }
    None
}

fn fasta_record<'a>(input: &'a [u8], last: bool) -> IResult<&[u8], SeqRecord<'a>> {
    let mut pos = 0;
    let id;
    match memchr(b'\n', &input) {
        None => return IResult::Incomplete(Needed::Unknown),
        Some(pos_end) => match str::from_utf8(&input[pos + 1..pos_end]) {
            Ok(v) => {
                id = v;
                pos += pos_end + 1;
            },
            // "Invalid UTF-8 in FASTA id: {}"
            Err(_) => return IResult::Error(ErrorKind::IsNotStr),
        },
    };

    match memchr2(b'\n', b'>', &input[pos..]) {
        None => match last {
            false => IResult::Incomplete(Needed::Unknown),
            true => IResult::Done(&[], (id, strip_whitespace(&input[pos..input.len()]), None)),
        },
        Some(pos_end) => {
            IResult::Done(&input[pos + pos_end + 1..],
                          (id, strip_whitespace(&input[pos..pos + pos_end]), None))
        },
    }
}

fn fastq_record<'a>(input: &'a [u8]) -> IResult<&[u8], SeqRecord<'a>> {
    if input[0] != b'@' {
        // TODO: return IResult::Error();
        panic!("Invalid FASTQ record");  // before line ...?
    }
    let mut pos = 0;
    let id;
    match memchr(b'\n', &input) {
        None => return IResult::Incomplete(Needed::Unknown),
        Some(pos_end) => match str::from_utf8(&input[pos + 1..pos_end]) {
            Ok(v) => {
                id = v;
                pos += pos_end + 1;
            },
            // "Invalid UTF-8 in FASTA id: {}"
            Err(_) => return IResult::Error(ErrorKind::IsNotStr),
        },
    };

    let seq;
    match memchr2(b'\n', b'+', &input[pos..]) {
        None => return IResult::Incomplete(Needed::Unknown),
        Some(pos_end) => {
            seq = input[pos..pos + pos_end].to_vec();
            pos += pos_end + 1;
        },
    };

    match memchr(b'\n', &input[pos..]) {
        None => return IResult::Incomplete(Needed::Unknown),
        Some(pos_end) => pos += pos_end + 1,
    };


    let end_len = pos + seq.len();
    if end_len > input.len() {
        IResult::Incomplete(Needed::Size(end_len - input.len()))
    } else {
        IResult::Done(&input[end_len + 1..], (id, seq, Some(&input[pos..end_len])))
    }
}

#[test]
fn test_parsers() {
    let parsed = fasta_record(b">\n", true);
    let res = ("", b"".to_vec(), None);
    assert_eq!(parsed, IResult::Done(&b""[..], res));

    let no_more = &b">"[..];

    let parsed = fasta_record(b">\n\n>", false);
    let res = ("", b"".to_vec(), None);
    assert_eq!(parsed, IResult::Done(no_more, res));

    let parsed = fasta_record(b">test\nagct\n>", false);
    let res = ("test", b"agct".to_vec(), None);
    assert_eq!(parsed, IResult::Done(no_more, res));

    let parsed = fasta_record(b">test2\nagct", true);
    let res = ("test2", b"agct".to_vec(), None);
    assert_eq!(parsed, IResult::Done(&b""[..], res));

    let parsed = fastq_record(b"@test\nagct\n+test\nAAAA\n");
    let res = ("test", b"agct".to_vec(), Some(&b"AAAA"[..]));
    assert_eq!(parsed, IResult::Done(&b""[..], res));
}


#[derive(PartialEq, Eq, Debug)]
enum FASTXState {
    Start,
    FASTA,
    FASTQ,
    Done,
}

struct FASTXConsumer<'x> {
    consumer_state: ConsumerState<(), (), Move>,
    state: FASTXState,
    callback: &'x mut for<'a> FnMut(SeqRecord<'a>) -> (),
}

impl<'x> FASTXConsumer<'x> {
    pub fn new<F>(callback: &'x mut F) -> FASTXConsumer<'x>
        where F: for<'a> FnMut(SeqRecord<'a>) -> (),
    {
        FASTXConsumer {
            consumer_state: ConsumerState::Continue(Move::Consume(0)),
            state: FASTXState::Start,
            callback: callback,
        }
    }
}

impl<'a, 'x> Consumer<&'a [u8], (), (), Move> for FASTXConsumer<'x> {
    fn state(&self) -> &ConsumerState<(), (), Move> {
        &self.consumer_state
    }

    fn handle(&mut self, input: Input<&[u8]>) -> &ConsumerState<(), (), Move> {
        match (&self.state, input) {
            (&FASTXState::Done, _) => {
                // TODO: read after close error
                self.consumer_state = ConsumerState::Error(());
            },
            (_, Input::Empty) |
            (_, Input::Eof(None)) => {
                // TODO: empty file error
                self.consumer_state = ConsumerState::Error(());
                self.state = FASTXState::Done;
            },
            (_, Input::Error) => {
                panic!("Error reading records; buffer may be too small");
            },
            (&FASTXState::Start, Input::Element(slice)) |
            (&FASTXState::Start, Input::Eof(Some(slice))) => {
                match slice[0] as char {
                    '>' => {
                        self.consumer_state = ConsumerState::Continue(Move::Consume(0usize));
                        self.state = FASTXState::FASTA;
                    },
                    '@' => {
                        self.consumer_state = ConsumerState::Continue(Move::Consume(0usize));
                        self.state = FASTXState::FASTQ;
                    },
                    _ => {
                        // TODO: not a valid FASTX file
                        self.consumer_state = ConsumerState::Error(());
                        self.state = FASTXState::Done;
                    },
                }
            },
            (&FASTXState::FASTA, inp) => {
                let (slice, last_record) = match inp {
                    Input::Element(slice) => (slice, false),
                    Input::Eof(Some(slice)) => (slice, true),
                    _ => panic!("Should never happen"),
                };
                match fasta_record(slice, last_record) {
                    IResult::Error(_) => {
                        // TODO: error report
                        self.consumer_state = ConsumerState::Error(());
                        self.state = FASTXState::Done;
                    },
                    IResult::Incomplete(n) => {
                        self.consumer_state = ConsumerState::Continue(Move::Await(n));
                    },
                    IResult::Done(remaining_slice, seq) => {
                        (self.callback)(seq);
                        if remaining_slice.len() == 0 {
                            self.state = FASTXState::Done;
                        }
                        self.consumer_state =
                            ConsumerState::Continue(Move::Consume(slice.len() -
                                                                  remaining_slice.len()));
                    },
                }
            },
            (&FASTXState::FASTQ, Input::Element(slice)) |
            (&FASTXState::FASTQ, Input::Eof(Some(slice))) => {
                match fastq_record(slice) {
                    IResult::Error(_) => {
                        // TODO: error report
                        self.consumer_state = ConsumerState::Error(());
                        self.state = FASTXState::Done;
                    },
                    IResult::Incomplete(n) => {
                        self.consumer_state = ConsumerState::Continue(Move::Await(n));
                    },
                    IResult::Done(remaining_slice, seq) => {
                        (self.callback)(seq);
                        if remaining_slice.len() <= 1 {
                            self.state = FASTXState::Done;
                        }
                        self.consumer_state =
                            ConsumerState::Continue(Move::Consume(slice.len() -
                                                                  remaining_slice.len()));
                    },
                }
            },
        }
        &self.consumer_state
    }
}


pub fn fastx_bytes<'b, F>(bytes: &'b [u8], ref mut callback: F)
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
{
    //! Parse a collection of bytes into FASTX records and calls `callback` on each.
    //!
    let mut producer = MemProducer::new(bytes, 10000000);
    let mut consumer = FASTXConsumer::new(callback);

    while consumer.state != FASTXState::Done {
        producer.apply(&mut consumer);
    }
}

pub fn fastx_file<F>(filename: &str, ref mut callback: F)
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
{
    //! Parse a file (given its name) into FASTX records and calls `callback` on each.
    //!
    //! Note there's currently a bug where the chunk size of 10Mb here prevents
    //! reading any individual FASTA records longer than that (and will cause the
    //! parser to prematurely end when it hits records longer than that).
    //! We should file a bug upstream on nom (or fix FileProducer's `fill` method to
    //! allocate longer than that).
    let mut producer = FileProducer::new(filename, 10000000).unwrap();
    let mut consumer = FASTXConsumer::new(callback);

    while consumer.state != FASTXState::Done {
        producer.apply(&mut consumer);
    }
}

#[test]
fn test_callback() {
    let mut i = 0;
    fastx_bytes(&b">test\nAGCT\n>test2\nGATC"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.0, "test");
                assert_eq!(&seq.1[..], &b"AGCT"[..]);
                assert_eq!(seq.2, None);
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(&seq.1[..], &b"GATC"[..]);
                assert_eq!(seq.2, None);
            },
            _ => {
                assert!(false);
            },
        }
        i += 1;
    });

    i = 0;
    fastx_file("./tests/data/test.fa", |seq| {
        match i {
            0 => {
                assert_eq!(seq.0, "test");
                assert_eq!(&seq.1[..], b"AGCTGATCGA");
                assert_eq!(seq.2, None);
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(&seq.1[..], b"TAGC");
                assert_eq!(seq.2, None);
            },
            _ => {
                assert!(false);
            },
        }
        i += 1;
    });
}

fn quality_mask<'a>(seq_rec: SeqRecord<'a>, ref score: u8) -> Vec<u8> {
    match seq_rec.2 {
        None => seq_rec.1,
        Some(quality) => {
            // could maybe speed this up by doing a copy of base and then
            // iterating though qual and masking?
            seq_rec.1
                   .iter()
                   .zip(quality.iter())
                   .map(|(base, qual)| {
                       if qual < score {
                           b'N'
                       } else {
                           base.clone()
                       }
                   })
                   .collect()

        },
    }
}

#[test]
fn test_quality_mask() {
    let seq_rec = ("", b"AGCT".to_vec(), Some(&b"AAA0"[..]));
    let filtered = quality_mask(seq_rec, '5' as u8);
    assert_eq!(&filtered[..], &b"AGCN"[..]);
}


#[test]
fn test_wrapped_fasta() {
    let mut i = 0;
    fastx_bytes(&b">test\nAGCT\nTCG\n>test2\nG"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.0, "test");
                assert_eq!(&seq.1[..], &b"AGCTTCG"[..]);
                assert_eq!(seq.2, None);
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(&seq.1[..], &b"G"[..]);
                assert_eq!(seq.2, None);
            },
            _ => {
                assert!(false);
            },
        }
        i += 1;
    });
}

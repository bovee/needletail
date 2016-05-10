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

use std::str;
use nom::{Err, ErrorKind, Input, IResult, Needed, Move};
use nom::{FileProducer, MemProducer, Consumer, ConsumerState, Producer};

/// A generic FASTX record containing:
///   seq.0 - an id
///   seq.1 - the sequence itself
///   seq.2 - an optional set of quality scores for the sequence
pub type SeqRecord<'a> = (&'a str, &'a [u8], Option<&'a [u8]>);

fn find_pos(buffer: &[u8], search: &[u8]) -> Option<usize> {
    //! Simple linear search through `buffer` looking for `search`
    // TODO: there has to be something like this built-in? can't find it though :/
    // TODO: (and in theory lots of room for optimization of longer
    // string cases here; but we should be okay for len 1)
    for (i, byte) in buffer.windows(search.len()).enumerate() {
        if *byte == *search {
            return Some(i);
        }
    }
    None
}

fn fasta_record<'a>(input: &'a [u8], last: bool) -> IResult<&[u8], SeqRecord<'a>> {
    let mut pos = 0;
    let id;
    match find_pos(&input, &['\n' as u8]) {
        None => return IResult::Incomplete(Needed::Unknown),
        Some(pos_end) => match str::from_utf8(&input[pos + 1..pos_end]) {
            Ok(v) => {
                id = v;
                pos += pos_end + 1;
            },
            Err(_) => return IResult::Error(Err::Code(ErrorKind::IsNotStr)), //"Invalid UTF-8 in FASTA id: {}"
        },
    };

    match find_pos(&input[pos..], &['\n' as u8, '>' as u8]) {
        None => match last {
            false => IResult::Incomplete(Needed::Unknown),
            true => IResult::Done(&[], (id, &input[pos..input.len()], None)),
        },
        Some(pos_end) => IResult::Done(&input[pos + pos_end + 1..], (id, &input[pos..pos + pos_end], None)),
    }
}

fn fastq_record<'a>(input: &'a [u8]) -> IResult<&[u8], SeqRecord<'a>> {
    let mut pos = 0;
    let id;
    match find_pos(&input, &['\n' as u8]) {
        None => return IResult::Incomplete(Needed::Unknown),
        Some(pos_end) => match str::from_utf8(&input[pos + 1..pos_end]) {
            Ok(v) => {
                id = v;
                pos += pos_end + 1;
            },
            Err(_) => return IResult::Error(Err::Code(ErrorKind::IsNotStr)), //"Invalid UTF-8 in FASTA id: {}"
        },
    };

    let seq;
    match find_pos(&input[pos..], &['\n' as u8, '+' as u8]) {
        None => return IResult::Incomplete(Needed::Unknown),
        Some(pos_end) => {
            seq = &input[pos..pos + pos_end];
            pos += pos_end + 1;
        },
    };

    match find_pos(&input[pos..], &['\n' as u8]) {
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
    let res = ("", &b""[..], None);
    assert_eq!(parsed, IResult::Done(&b""[..], res));
    
    let no_more = &b">"[..];

    let parsed = fasta_record(b">\n\n>", false);
    let res = ("", &b""[..], None);
    assert_eq!(parsed, IResult::Done(no_more, res));

    let parsed = fasta_record(b">test\nagct\n>", false);
    let res = ("test", &b"agct"[..], None);
    assert_eq!(parsed, IResult::Done(no_more, res));

    let parsed = fasta_record(b">test2\nagct", true);
    let res = ("test2", &b"agct"[..], None);
    assert_eq!(parsed, IResult::Done(&b""[..], res));

    let parsed = fastq_record(b"@test\nagct\n+test\nAAAA\n");
    let res = ("test", &b"agct"[..], Some(&b"AAAA"[..]));
    assert_eq!(parsed, IResult::Done(&b""[..], res));
}


#[derive(PartialEq, Eq, Debug)]
enum FASTXState {
    Start,
    FASTA,
    FASTQ,
    Done,
}

struct FASTXConsumer<'x>  {
    consumer_state: ConsumerState<(), (), Move>,
    state: FASTXState,
    callback: &'x mut for<'a> FnMut(SeqRecord<'a>) -> (),
}

impl<'x> FASTXConsumer<'x> {
    pub fn new<F>(callback: &'x mut F) -> FASTXConsumer<'x> where F: for<'a> FnMut(SeqRecord<'a>) -> ()  {
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
            (_, Input::Empty) | (_, Input::Eof(None)) => {
                // TODO: empty file error
                self.consumer_state = ConsumerState::Error(());
                self.state = FASTXState::Done;
            },
            (&FASTXState::Start, Input::Element(sl)) | (&FASTXState::Start, Input::Eof(Some(sl))) => {
                match sl[0] as char {
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
                let (sl, last_record) = match inp {
                    Input::Element(sl) => (sl, false),
                    Input::Eof(Some(sl)) => (sl, true),
                    _ => panic!("Should never happen"),
                };
                match fasta_record(sl, last_record) {
                    IResult::Error(_) => {
                        // TODO: error report
                        self.consumer_state = ConsumerState::Error(());
                        self.state = FASTXState::Done;
                    },
                    IResult::Incomplete(n) => {
                        self.consumer_state = ConsumerState::Continue(Move::Await(n));
                    },
                    IResult::Done(remaining_sl, seq) => {
                        (self.callback)(seq);
                        if remaining_sl.len() == 0 {
                            self.state = FASTXState::Done;
                        }
                        self.consumer_state =  ConsumerState::Continue(Move::Consume(sl.len() - remaining_sl.len()));
                    }
                }
            },
            (&FASTXState::FASTQ, Input::Element(sl)) | (&FASTXState::FASTQ, Input::Eof(Some(sl))) => {
                match fastq_record(sl) {
                    IResult::Error(_) => {
                        // TODO: error report
                        self.consumer_state = ConsumerState::Error(());
                        self.state = FASTXState::Done;
                    },
                    IResult::Incomplete(n) => {
                        self.consumer_state = ConsumerState::Continue(Move::Await(n));
                    },
                    IResult::Done(remaining_sl, seq) => {
                        (self.callback)(seq);
                        if remaining_sl.len() <= 1 {
                            self.state = FASTXState::Done;
                        }
                        self.consumer_state =  ConsumerState::Continue(Move::Consume(sl.len() - remaining_sl.len()));
                    },
                }
            },
        }
        &self.consumer_state
    }
}


pub fn fastx_bytes<'b, F>(bytes: &'b [u8], ref mut callback: F) where F: for<'a> FnMut(SeqRecord<'a>) -> () {
    //! Parse a collection of bytes into FASTX records and calls `callback` on each.
    //!
    let mut producer = MemProducer::new(bytes, 10000000);
    let mut consumer = FASTXConsumer::new(callback);

    while consumer.state != FASTXState::Done {
        producer.apply(&mut consumer);
    }
}

pub fn fastx_file<F>(filename: &str, ref mut callback: F) where F: for<'a> FnMut(SeqRecord<'a>) -> () {
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
                assert_eq!(seq.1, &b"AGCT"[..]);
                assert_eq!(seq.2, None);
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(seq.1, &b"GATC"[..]);
                assert_eq!(seq.2, None);
            },
            _ => {
                assert!(false);
            }
        }
        i += 1;
    });

    i = 0;
    fastx_file("./test.fa", |seq| {
        match i {
            0 => {
                assert_eq!(seq.0, "test");
                assert_eq!(seq.1, &b"AGCTGATCGA"[..]);
                assert_eq!(seq.2, None);
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(seq.1, &b"TAGC\n"[..]);
                assert_eq!(seq.2, None);
            },
            _ => {
                assert!(false);
            }
        }
        i += 1;
    });
}
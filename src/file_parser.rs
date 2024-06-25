pub mod fileparser {
    use std::char;
    use std::io::Error;
    use std::path::Path;

    pub struct DNASequence {
        pub description: String,
        pub record: String,
        pub iscircular: bool,
        pub features: Vec<Feature>,
        pub iter_pos: usize,
        pub max_iter_size: usize,
        pub buffered_sequence: String,
        pub sequence_length: usize,
        buffer_is_full: bool,
    }

    impl DNASequence {
        fn new(seq: String, max_iter_size: usize, iscircular: bool) -> Result<DNASequence, Error> {
            // Perform sanity checking

            if max_iter_size <= 0 {
                return Err(Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "The max_iter_size must be greater than 0",
                ));
            }

            if max_iter_size > seq.chars().count() {
                return Err(Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "The max_iter_size must be less than or equal to the length of the sequence",
                ));
            }

            let valid_bases = "ACGTURYSWKMBDHVNacgturyswkmbdhvn";
            for base in seq.chars() {
                if !valid_bases.contains(base) {
                    return Err(Error::new(
                        std::io::ErrorKind::InvalidInput,
                        "Invalid base in sequence",
                    ));
                }
            }

            // Make the new sequence and return it
            let result: DNASequence;
            unsafe { result = Self::new_unchecked(seq, max_iter_size, iscircular) }
            Ok(result)
        }

        /// Creates a new 'DNAsequence' for parsing DNA. Non-valid base characters in seq may lead to nondeterministic behaviour
    unsafe fn new_unchecked(
        seq: String,
        max_iter_size: usize,
        iscircular: bool,
    ) -> DNASequence {
        let seq_length = &seq.chars().count();
        DNASequence {
            description: String::new(),
            record: seq,
            iscircular: iscircular,
            features: Vec::new(),
            iter_pos: 0,
            max_iter_size: max_iter_size,
            buffered_sequence: String::new(),
            sequence_length: seq_length.clone(),
            buffer_is_full: false,
        }
     }
    }

    impl Iterator for DNASequence {
        type Item = SeqSegment;

        fn next(&mut self) -> Option<SeqSegment> {
            // If we're at the end of the sequence, return none
            if self.iter_pos >= self.sequence_length {
                return None;
            } // TODO: Loop around for circular DNA

            // Add the next character to the buffered sequence until the max length is reached

            if self.buffer_is_full {
                // If the buffered sequence is full, remove the first character and add the next
                self.buffered_sequence.remove(0);
                self.buffered_sequence
                    .push(self.record.as_bytes()[self.iter_pos as usize] as char);
                self.iter_pos += 1;
                Some(SeqSegment {
                    sequence: self.buffered_sequence.clone(),
                    start: self.iter_pos - self.max_iter_size,
                    end: self.iter_pos,
                })
            } else if !self.buffer_is_full
                && self.buffered_sequence.chars().count() < self.max_iter_size
            {
                // If buffer isn't full, add the next
                self.buffered_sequence
                    .push(self.record.as_bytes()[self.iter_pos as usize] as char);
                self.iter_pos += 1;
                Some(SeqSegment {
                    sequence: self.buffered_sequence.clone(),
                    start: self.iter_pos - self.buffered_sequence.chars().count(),
                    end: self.iter_pos,
                })
            } else {
                // If this fills the buffer, remove the first character and add the next, then set the buffer to full
                self.buffered_sequence.remove(0);
                self.buffered_sequence
                    .push(self.record.as_bytes()[self.iter_pos as usize] as char);
                self.iter_pos += 1;
                self.buffer_is_full = true;
                Some(SeqSegment {
                    sequence: self.buffered_sequence.clone(),
                    start: self.iter_pos - self.max_iter_size,
                    end: self.iter_pos,
                })
            }
        }
    }

    pub struct Feature {
        description: String,
        start: usize,
        end: usize,
    }

    fn parse_fasta(
        file: &Path,
        max_iter_size: usize,
        iscircular: bool,
    ) -> Result<Vec<DNASequence>, std::io::Error> {
        use std::fs::File;
        use std::io::prelude::*;
        use std::io::BufReader;
        let file = File::open(file)?;
        let mut f = BufReader::new(file);
        let mut buf = Vec::<u8>::new();
        let mut description = String::new();
        let mut sequence = String::new();

        let mut sequences = Vec::new();

        let mut add_to_description = false;
        let mut in_seq = false;

        // TODO: Add support for multifasta
        while f.read_until(b'\n', &mut buf).expect("read_until failed") != 0 {
            let s = String::from_utf8(buf).expect("from_utf8 failed");
            for c in s.chars() {
                if c == '>' {
                    if in_seq {
                        in_seq = false;
                        let new_seq = DNASequence::new(sequence, max_iter_size, iscircular);
                        match new_seq {
                            Ok(seq) => sequences.push(seq),
                            Err(e) => return Err(e),
                        };
                        description = String::new();
                        sequence = String::new();
                        add_to_description = true;
                    } else {
                        add_to_description = true;
                    }
                } else if c == '\n' && !in_seq {
                    add_to_description = false;
                    in_seq = true;
                } else if add_to_description {
                    description.push(c);
                } else if c.is_alphabetic() && in_seq {
                    sequence.push(c);
                }
            }
            // this returns the ownership of the read data to buf
            // there is no allocation
            buf = s.into_bytes();
            buf.clear();

            if !in_seq && !add_to_description {
                break;
            }
        }

        if in_seq {
            let new_seq = DNASequence::new(sequence, max_iter_size, iscircular);
            match new_seq {
                Ok(seq) => sequences.push(seq),
                Err(e) => return Err(e),
            };
        }

        Ok(sequences)
    }

    struct GenbankParser {
        name: String,
        file: String,
        iscircular: bool,
        bufferedseq: String,
    }
    impl GenbankParser {
        fn new(_file: String, _iscircular: bool) {
            todo!("Implement GenbankParser")
        }

        fn get(self, _start: i32, _end: i32) -> SeqSegment {
            todo!("Implement get")
        }
        fn next() -> SeqSegment {
            todo!("Implement next")
        }
    }

    pub struct SeqSegment {
        pub sequence: String,
        pub start: usize,
        pub end: usize,
    }

    pub fn parse_file(
        filename: &str,
        max_iter_size: usize,
    ) -> Result<Vec<DNASequence>, std::io::Error> {
        // Check the file type
        // Create the appropriate SeqParser
        // Return the SeqParser

        let iscircular = true;
        let file = Path::new(filename);

        // Check to see if the file exists
        if !file.exists() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "File does not exist",
            ));
        }

        println!("Parsing file");

        if file.extension().unwrap() == "fasta" {
            let result = parse_fasta(file, max_iter_size, iscircular);
            match result {
                Ok(seq) => return Ok(seq),
                Err(e) => return Err(e),
            }
        } else {
            // Raise an error
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "File type not supported",
            ))
        }
    }
}

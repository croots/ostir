use std::path::Path;

/// Parse multi-FASTA text, returning (name, sequence) pairs.
pub fn parse_fasta_str(content: &str) -> Vec<(String, String)> {
    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = String::new();

    for line in content.lines() {
        let line = line.trim();
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                sequences.push((current_name.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_name = line[1..].trim().to_string();
            if current_name.is_empty() {
                current_name = "unnamed".to_string();
            }
        } else if !line.is_empty() {
            current_seq.push_str(line);
        }
    }
    if !current_seq.is_empty() {
        sequences.push((current_name, current_seq));
    }
    sequences
}

/// A single input sequence with all parameters for an OSTIR run.
#[derive(Debug, Clone)]
pub struct OstirInput {
    pub sequence: String,
    pub name: String,
    pub asd: String,
    pub start: usize,  // 1-indexed
    pub end: Option<usize>, // 1-indexed; None = sequence length
    pub circular: bool,
    pub print_sequence: bool,
    pub print_asd: bool,
}

impl OstirInput {
    pub fn new(sequence: String, name: String) -> Self {
        OstirInput {
            sequence,
            name,
            asd: "ACCTCCTTA".to_string(),
            start: 1,
            end: None,
            circular: false,
            print_sequence: false,
            print_asd: false,
        }
    }
}

/// Parse a FASTA file into OstirInput records, applying CLI defaults.
pub fn parse_fasta_file(path: &Path, defaults: &OstirInput) -> Result<Vec<OstirInput>, std::io::Error> {
    let content = std::fs::read_to_string(path)?;
    let pairs = parse_fasta_str(&content);
    Ok(pairs
        .into_iter()
        .map(|(name, seq)| OstirInput {
            sequence: seq,
            name,
            asd: defaults.asd.clone(),
            start: defaults.start,
            end: defaults.end,
            circular: defaults.circular,
            print_sequence: defaults.print_sequence,
            print_asd: defaults.print_asd,
        })
        .collect())
}

/// Parse a CSV file into OstirInput records.
/// Supports columns: sequence/seq, name/id, anti-Shine-Dalgarno, start, end, circular
pub fn parse_csv_file(path: &Path, defaults: &OstirInput) -> Result<Vec<OstirInput>, std::io::Error> {
    let content = std::fs::read_to_string(path)?;
    let mut inputs = Vec::new();
    let mut idx = 0usize;

    // Strip blank/comment lines then feed to csv reader
    let filtered: String = content
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.trim().starts_with('#'))
        .map(|l| format!("{}\n", l))
        .collect();

    let mut rdr = csv::Reader::from_reader(filtered.as_bytes());
    let headers: Vec<String> = rdr
        .headers()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?
        .iter()
        .map(|h| h.to_lowercase())
        .collect();

    for record in rdr.records() {
        idx += 1;
        let record = record.map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
        let get = |key: &str| -> Option<String> {
            headers.iter().position(|h| h == key).and_then(|i| record.get(i)).map(|v| v.to_string()).filter(|v| !v.is_empty())
        };

        let sequence = get("sequence").or_else(|| get("seq"))
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing 'sequence'/'seq' column"))?;
        let sequence = sequence.replace(' ', "");

        let name = get("name").or_else(|| get("id")).unwrap_or_else(|| format!("sequence_{}", idx));

        let asd = get("anti-shine-dalgarno")
            .or_else(|| get("anti_shine_dalgarno"))
            .unwrap_or_else(|| defaults.asd.clone());
        let asd = if asd.is_empty() { defaults.asd.clone() } else { asd };

        let start = get("start")
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(defaults.start);
        let end = get("end")
            .and_then(|v| v.parse::<usize>().ok())
            .or(defaults.end);
        let circular = get("circular")
            .map(|v| matches!(v.to_lowercase().as_str(), "true" | "t" | "yes" | "y" | "1"))
            .unwrap_or(defaults.circular);

        inputs.push(OstirInput {
            sequence,
            name,
            asd,
            start,
            end,
            circular,
            print_sequence: defaults.print_sequence,
            print_asd: defaults.print_asd,
        });
    }
    Ok(inputs)
}

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

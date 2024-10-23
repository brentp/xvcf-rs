use rust_htslib::bcf::{self, Read};
use std::{fmt, io, path::Path};

pub enum ReaderInner {
    Indexed(bcf::IndexedReader),
    Plain(bcf::Reader),
}

pub struct Reader<'a> {
    inner: ReaderInner,
    header: bcf::header::HeaderView,
    records: Option<bcf::Records<'a, bcf::Reader>>,
    current_record: Option<bcf::Record>,
}

impl<'a> Reader<'a> {
    pub fn new(inner: ReaderInner) -> Self {
        let header = match &inner {
            ReaderInner::Indexed(r) => r.header().clone(),
            ReaderInner::Plain(r) => r.header().clone(),
        };
        Self {
            inner,
            header,
            records: None,
            current_record: None,
        }
    }

    pub fn take(&mut self) -> Option<bcf::Record> {
        self.current_record.take()
    }

    pub fn next_record(&mut self) -> io::Result<Option<bcf::Record>> {
        if let Some(record) = self.current_record.take() {
            return Ok(Some(record));
        }

        match &mut self.inner {
            ReaderInner::Indexed(reader) => {
                let mut record = reader.empty_record();
                if let Some(records) = self.records.as_mut() {
                    return records
                        .next()
                        .transpose()
                        .map_err(|e| io::Error::new(io::ErrorKind::Other, e));
                }
                match reader.read(&mut record) {
                    Some(Ok(())) => Ok(Some(record)),
                    Some(Err(e)) => Err(io::Error::new(io::ErrorKind::Other, e)),
                    None => Ok(None),
                }
            }
            ReaderInner::Plain(reader) => match reader.records().next() {
                Some(Ok(record)) => Ok(Some(record)),
                Some(Err(e)) => Err(io::Error::new(io::ErrorKind::Other, e)),
                None => Ok(None),
            },
        }
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Reader<'static>> {
        // Try indexed reader first
        match bcf::IndexedReader::from_path(&path) {
            Ok(indexed_reader) => Ok(Reader::new(ReaderInner::Indexed(indexed_reader))),
            Err(_) => {
                // Fall back to plain reader
                let reader = bcf::Reader::from_path(path)
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
                Ok(Reader::new(ReaderInner::Plain(reader)))
            }
        }
    }

    pub fn header(&self) -> &bcf::header::HeaderView {
        &self.header
    }
}

pub trait Skip {
    fn skip_to(&mut self, region: &str) -> io::Result<()>;
}

impl<'a> Skip for Reader<'a> {
    fn skip_to(&mut self, region: &str) -> io::Result<()> {
        // Parse region string (e.g., "chr1:2000")
        let parts: Vec<&str> = region.split(':').collect();
        if parts.len() != 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid region format",
            ));
        }

        let chrom = parts[0];
        let pos: u64 = parts[1]
            .parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        match &mut self.inner {
            ReaderInner::Indexed(reader) => {
                let rid = reader
                    .header()
                    .name2rid(chrom.as_bytes())
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

                reader
                    .fetch(rid, pos - 1, None)
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))
            }
            ReaderInner::Plain(reader) => {
                let target_rid = reader
                    .header()
                    .name2rid(chrom.as_bytes())
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

                // Scan through records until we find one that's >= our target position
                match reader.records().next() {
                    Some(Ok(record)) => {
                        let record_rid = record.rid().unwrap();

                        if record_rid > target_rid {
                            // We've gone past our target chromosome
                            self.current_record = Some(record);
                            return Ok(());
                        } else if record_rid == target_rid {
                            if (record.end() as u64) >= pos {
                                // Found a record at or after our target position
                                self.current_record = Some(record);
                                return Ok(());
                            }
                        }
                        // Otherwise continue scanning
                    }
                    Some(Err(e)) => return Err(io::Error::new(io::ErrorKind::Other, e)),
                    None => return Ok(()),
                }
                Ok(()) // Reached end of file without finding matching position
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_read_vcf() {
        let path = PathBuf::from("tests/t.vcf.gz");
        let mut reader = Reader::from_path(&path).expect("error creating new reader");
        let _header = reader.header().clone();

        reader
            .skip_to("chr1:2000")
            .expect("error skipping to region");

        // Store the record result in a separate variable to end the borrow
        let record = reader.next_record().expect("error reading record");
        if let Some(record) = record {
            assert_eq!(record.pos(), 1999); // 0-based position
            let chrom = reader.header().rid2name(record.rid().unwrap()).unwrap();
            let chrom_str = std::str::from_utf8(chrom).unwrap();
            println!("Record: {:?}", record);
            println!("Record: {}:{}", chrom_str, record.pos());
        } else {
            panic!("No record found");
        }

        // now if we skip backwards we should NOT get the same record
        reader
            .skip_to("chr1:1999")
            .expect("error skipping to region");
        let r = reader.next_record().expect("error reading record").unwrap();
        let chrom = reader.header().rid2name(r.rid().unwrap()).unwrap();
        let chrom_str = std::str::from_utf8(chrom).unwrap();
        eprintln!("Record: {}:{}", chrom_str, r.pos());
        assert!(r.pos() != 1999);
    }
}

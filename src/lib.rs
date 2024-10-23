use rust_htslib::bcf::{self, Read};
use std::{io, path::Path};

pub enum ReaderInner {
    Indexed(bcf::IndexedReader),
    Plain(bcf::Reader),
}

pub struct Reader<'a> {
    inner: ReaderInner,
    header: bcf::header::HeaderView,
    records: Option<bcf::Records<'a, bcf::Reader>>,
    current_record: Option<bcf::Record>,
    last_record: TinyRecord,
}

#[derive(Clone, Debug)]
struct TinyRecord {
    rid: i32,
    pos: i64,
    stop: i64,
}

impl Reader<'_> {
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
            last_record: TinyRecord {
                rid: -1,
                pos: -1,
                stop: -1,
            },
        }
    }

    pub fn take(&mut self) -> Option<bcf::Record> {
        if let Some(record) = self.current_record.take() {
            self.last_record = TinyRecord {
                rid: record.rid().unwrap() as i32,
                pos: record.pos(),
                stop: record.end(),
            };
            Some(record)
        } else {
            None
        }
    }

    pub fn next_record(&mut self) -> io::Result<Option<bcf::Record>> {
        if let Some(record) = self.take() {
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
        let path = path.as_ref();

        // Check for .csi or .tbi index files
        let has_index = Path::new(&format!("{}.csi", path.display())).exists()
            || Path::new(&format!("{}.tbi", path.display())).exists();

        if has_index {
            // Use indexed reader if index exists
            match bcf::IndexedReader::from_path(path) {
                Ok(indexed_reader) => Ok(Reader::new(ReaderInner::Indexed(indexed_reader))),
                Err(e) => {
                    // Fall back to plain reader if index exists but can't be loaded
                    eprintln!(
                        "Index found but failed to load: {}. Falling back to plain reader",
                        e
                    );
                    let reader = bcf::Reader::from_path(path)
                        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
                    Ok(Reader::new(ReaderInner::Plain(reader)))
                }
            }
        } else {
            // Use plain reader if no index exists
            let reader = bcf::Reader::from_path(path)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            Ok(Reader::new(ReaderInner::Plain(reader)))
        }
    }

    pub fn header(&self) -> &bcf::header::HeaderView {
        &self.header
    }
}

fn is_record_after_last(last: &TinyRecord, record: &bcf::Record) -> bool {
    let current_rid = record.rid().unwrap();

    if current_rid as i32 > last.rid {
        return true;
    }
    if (current_rid as i32) < last.rid {
        return false;
    }
    if record.pos() < last.pos {
        return false;
    }
    if record.pos() > last.pos {
        return true;
    }
    record.end() > last.stop
}

pub trait Skip {
    fn skip_to(&mut self, region: &str) -> io::Result<()>;
}

impl Skip for Reader<'_> {
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
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

                // Read until we find a record after the last_record
                let mut record = reader.empty_record();
                loop {
                    match reader.read(&mut record) {
                        Some(Ok(())) => {
                            if is_record_after_last(&self.last_record, &record) {
                                /*
                                eprintln!(
                                    "setting current record to: {} with last_record: {:?}",
                                    record.pos(),
                                    self.last_record
                                );
                                */
                                self.current_record = Some(record);
                                return Ok(());
                            }
                        }
                        Some(Err(e)) => return Err(io::Error::new(io::ErrorKind::Other, e)),
                        None => {
                            break;
                        }
                    }
                }
                // if we got here, then we probably hit the end of a chrom. so skip to next chrom.
                let next_rid = rid + 1;
                let name = self.header().rid2name(next_rid);
                if let Err(_e) = name {
                    return Ok(());
                }
                let name: String = std::str::from_utf8(name.unwrap()).unwrap().to_owned() + ":1";
                self.skip_to(&name)
            }
            ReaderInner::Plain(reader) => {
                let target_rid = reader
                    .header()
                    .name2rid(chrom.as_bytes())
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

                // Scan through records until we find one that's >= our target position
                // AND after our last_record
                loop {
                    match reader.records().next() {
                        Some(Ok(record)) => {
                            let record_rid = record.rid().unwrap();

                            if record_rid > target_rid {
                                // We've gone past our target chromosome
                                if is_record_after_last(&self.last_record, &record) {
                                    self.current_record = Some(record);
                                }
                                return Ok(());
                            } else if record_rid == target_rid
                                && (record.end() as u64) >= pos
                                && is_record_after_last(&self.last_record, &record)
                            {
                                // Found a record at or after our target position
                                self.current_record = Some(record);
                                return Ok(());
                            }
                            // Otherwise continue scanning
                        }
                        Some(Err(e)) => return Err(io::Error::new(io::ErrorKind::Other, e)),
                        None => return Ok(()),
                    }
                }
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
            println!("first seek: {}:{}", chrom_str, record.pos());
        } else {
            panic!("No record found");
        }

        // now if we skip backwards we should NOT get the same record
        reader
            .skip_to("chr1:2000")
            .expect("error skipping to region");
        let record = reader.next_record().expect("error reading record");
        assert!(record.is_some());
        let record = record.unwrap();
        assert_eq!(record.pos(), 2999);
        assert_eq!(record.rid().unwrap(), 1);
    }
}

pub use noodles::bcf;
use noodles::bgzf;
use noodles::core::{Position, Region};
pub use noodles::csi;
use noodles::vcf::record::Chromosome;
use noodles::vcf::IndexedReader;
pub use noodles::vcf::{self};

pub use noodles_util::variant;

pub mod detect;
use detect::{Compression, Format};
use std::io::{self, Seek};
use std::io::{BufRead, BufReader};

pub enum XCF {
    /*
    Vcf(vcf::Reader<BufReader<Box<dyn BufRead>>>),
    Bcf(bcf::Reader<BufReader<Box<dyn BufRead>>>),
    IndexedVcf(vcf::indexed_reader::IndexedReader<BufReader<Box<dyn BufRead>>>),
    IndexedBcf(bcf::indexed_reader::IndexedReader<BufReader<Box<dyn BufRead>>>),
    */
    Vcf(vcf::Reader<BufReader<Box<dyn BufRead>>>),
    Bcf(bcf::Reader<BufReader<Box<dyn BufRead>>>),
    IndexedVcf(vcf::indexed_reader::IndexedReader<bgzf::Reader<BufReader<Box<dyn BufRead>>>>),
    IndexedBcf(bcf::indexed_reader::IndexedReader<bgzf::Reader<BufReader<Box<dyn BufRead>>>>),
    CompressedBcf(bcf::Reader<bgzf::Reader<BufReader<Box<dyn BufRead>>>>),
    CompressedVcf(vcf::Reader<bgzf::Reader<BufReader<Box<dyn BufRead>>>>),
}

pub struct Reader {
    inner: XCF,
    header: vcf::Header,
    variant: Option<vcf::Record>,
}

impl Reader {
    pub fn new(inner: XCF, header: vcf::Header) -> Self {
        Self {
            inner,
            header,
            variant: None,
        }
    }

    pub fn from_reader(reader: Box<dyn BufRead>, path: Option<String>) -> io::Result<Reader> {
        let mut reader = BufReader::new(reader);
        let compression = detect::detect_compression(&mut reader)?;
        let format = detect::detect_format(&mut reader, compression)?;
        let csi = find_index(path);

        let mut rdr = match (format, compression) {
            (Format::Vcf, None) => {
                let mut reader = vcf::Reader::new(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(reader), header)
            }
            (Format::Vcf, Some(Compression::Bgzf)) => {
                let mut bgzf_reader = noodles::bgzf::Reader::new(reader);
                if let Some(csi) = csi {
                    let mut reader = IndexedReader::new(bgzf_reader, csi);
                    let header = reader.read_header()?;
                    Reader::new(XCF::IndexedVcf(reader), header)
                } else {
                    let mut reader = vcf::Reader::new(bgzf_reader);
                    let header = reader.read_header()?;
                    Reader::new(XCF::CompressedVcf(reader), header)
                }
            }
            (Format::Bcf, None) => {
                let mut reader = bcf::Reader::from(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Bcf(reader), header)
            }
            (Format::Bcf, Some(Compression::Bgzf)) => {
                if let Some(csi) = csi {
                    let mut reader = bcf::IndexedReader::new(reader, csi);
                    let header = reader.read_header()?;
                    Reader::new(XCF::IndexedBcf(reader), header)
                } else {
                    let mut bgzf_reader = noodles::bgzf::Reader::new(reader);
                    let mut reader = bcf::Reader::from(bgzf_reader);
                    let header = reader.read_header()?;
                    Reader::new(XCF::CompressedBcf(reader), header)
                }
            }
        };
        Ok(rdr)
    }

    pub fn read_record(&mut self, header: &vcf::Header, v: &mut vcf::Record) -> io::Result<usize> {
        if let Some(variant) = self.variant.take() {
            *v = variant;
            return Ok(1);
        }
        match &mut self.inner {
            XCF::Vcf(reader) => reader.read_record(header, v),
            XCF::Bcf(reader) => reader.read_record(header, v),
            XCF::IndexedVcf(reader) => reader.read_record(header, v),
            XCF::IndexedBcf(reader) => reader.read_record(header, v),
            XCF::CompressedBcf(reader) => reader.read_record(header, v),
            XCF::CompressedVcf(reader) => reader.read_record(header, v),
        }
    }

    pub fn header(&mut self) -> &vcf::Header {
        &self.header
    }
}

/*
#[inline]
fn chrom_equals(c: &Chromosome, name: &str) -> bool {
    match c {
        Chromosome::Name(n) => n == name,
        Chromosome::Symbol(_) => false,
    }
}

fn advance_reader(r: &mut Reader, header: &vcf::Header, region: &Region) -> io::Result<()>
where
    R: BufRead,
{
    let one = Position::try_from(1).unwrap();
    let start =
        vcf::record::Position::try_from(usize::from(region.interval().start().unwrap_or(one)))
            .unwrap();
    let mut v = vcf::Record::default();
    r.variant = None;
    loop {
        let result = match &mut r.inner {
            XCF::Vcf(reader) => reader.read_record(header, &mut v),
            XCF::Bcf(reader) => reader.read_record(header, &mut v),
            _ => unimplemented!(),
        };
        match result {
            Ok(0) => return Ok(()),
            Ok(_) => {
                if chrom_equals(v.chromosome(), region.name()) && v.position() >= start {
                    r.variant = Some(v);
                    return Ok(());
                }
            }
            Err(e) => return Err(e),
        };
    }
}

impl<R> Reader<R>
where
    R: BufRead + Seek,
{
    // skip_to simply sets the file pointer to the start of region.
    // internally, it consumes the first variant, but that will be returned on the
    // first call to read_record.
    pub fn skip_to(&mut self, header: &vcf::Header, region: &Region) -> io::Result<()> {
        match &mut self.inner {
            XCF::Vcf(_) => advance_reader(self, header, region),
            XCF::Bcf(_) => advance_reader(self, header, region),
            XCF::IndexedVcf(reader) => match reader.query(header, region) {
                Ok(mut r) => {
                    self.variant = match r.next() {
                        Some(Ok(v)) => Some(v),
                        Some(Err(e)) => return Err(e),
                        None => None,
                    };
                    Ok(())
                }
                Err(e) => Err(e),
            },
            XCF::IndexedBcf(reader) => match reader.query(header, region) {
                Ok(mut r) => {
                    self.variant = match r.next() {
                        Some(Ok(v)) => Some(v),
                        Some(Err(e)) => return Err(e),
                        None => None,
                    };
                    Ok(())
                }
                Err(e) => Err(e),
            },
        }
    }
}
*/

fn find_index(path: Option<String>) -> Option<csi::Index> {
    if let Some(path) = path {
        let mut index_path = path.clone();
        // TODO: look for csi or .tbi or use ##idx## in path.
        index_path.push_str(".csi");
        csi::read(index_path)
            .or_else(|_| {
                let mut index_path = path.clone();
                index_path.push_str(".tbi");
                csi::read(index_path)
            })
            .ok()
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_read_vcf() {
        let vcf_data = b"\
            ##fileformat=VCFv4.2\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
            chr1\t1000\t.\tA\tC\t.\t.\t.\n\
            chr1\t2000\t.\tG\tT\t.\t.\t.\n\
            chr2\t3000\t.\tC\tG\t.\t.\t.\n\
        ";
        let cursor = Cursor::new(vcf_data);
        let mut rdr = BufReader::new(std::fs::File::open("tests/t.bcf").unwrap());
        rdr.fill_buf().unwrap();
        eprintln!("buffer: {:?}", rdr.buffer());
        let rdr = Box::new(rdr);

        let mut reader = Reader::from_reader(rdr, None).expect("error creating new reader");
        let header = reader.header().clone();

        let start = Position::try_from(2000).expect("error creating start");
        let stop = Position::try_from(2100).expect("error creating stop");
        let region = Region::new("chr1", start..=stop);
        //reader.skip_to(&header, &region).unwrap();

        let mut v = vcf::Record::default();
        while let Ok(_) = reader.read_record(&header, &mut v) {
            eprintln!("v: {:?}", v);
            assert!(chrom_equals(v.chromosome(), "chr1"));
            assert_eq!(v.position(), start);
            break;
        }
    }
}

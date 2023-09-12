pub use noodles::bcf;
use noodles::core::{Position, Region};
pub use noodles::csi;
use noodles::vcf::record::Chromosome;
pub use noodles::vcf::{self};

pub use noodles_util::variant;

pub mod detect;
use detect::{Compression, Format};
use std::io::{self, Seek};
use std::io::{BufRead, BufReader};

pub enum XCF<R> {
    Vcf(vcf::Reader<R>),
    Bcf(bcf::Reader<R>),
    IndexedVcf(vcf::IndexedReader<R>),
    IndexedBcf(vcf::IndexedReader<R>),
}

pub struct Reader<R> {
    inner: XCF<R>,
    index: Option<csi::Index>,
    header: vcf::Header,
    variant: Option<vcf::Record>,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: XCF<R>, index: Option<csi::Index>, header: vcf::Header) -> Self {
        Self {
            inner,
            index,
            header,
            variant: None,
        }
    }
    pub fn from_reader(reader: R, path: Option<String>) -> io::Result<Reader<BufReader<R>>>
    where
        R: BufRead + 'static,
    {
        let mut reader = BufReader::new(reader);
        let compression = detect::detect_compression(&mut reader)?;
        let format = detect::detect_format(&mut reader, compression)?;

        let mut rdr = match (format, compression) {
            (Format::Vcf, None) => {
                let mut reader = vcf::Reader::new(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(reader), None, header)
            }
            (Format::Vcf, Some(Compression::Bgzf)) => {
                let mut reader = vcf::Reader::new(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(reader), None, header)
            }
            (Format::Bcf, None) => {
                let mut reader = bcf::Reader::from(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Bcf(reader), None, header)
            }
            (Format::Bcf, Some(Compression::Bgzf)) => {
                let mut reader = bcf::Reader::from(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Bcf(reader), None, header)
            }
        };
        rdr.index = find_index(path);
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
        }
    }

    pub fn header(&mut self) -> &vcf::Header {
        &self.header
    }

    pub fn index(&mut self) -> Option<&csi::Index> {
        self.index.as_ref()
    }
}

#[inline]
fn chrom_equals(c: &Chromosome, name: &str) -> bool {
    match c {
        Chromosome::Name(n) => n == name,
        Chromosome::Symbol(_) => false,
    }
}

fn advance_reader<R>(r: &mut Reader<R>, header: &vcf::Header, region: &Region) -> io::Result<()>
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
        let mut reader = Reader::from_reader(cursor, None).expect("error creating new reader");
        let header = reader.header().clone();

        let start = Position::try_from(2000).expect("error creating start");
        let stop = Position::try_from(2100).expect("error creating stop");
        let region = Region::new("chr1", start..=stop);
        reader.skip_to(&header, &region).unwrap();

        let mut v = vcf::Record::default();
        while let Ok(_) = reader.read_record(&header, &mut v) {
            eprintln!("v: {:?}", v);
            assert!(chrom_equals(v.chromosome(), "chr1"));
            assert_eq!(v.position(), start);
            break;
        }
    }
}

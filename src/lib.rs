pub use noodles::bcf;
use noodles::bgzf;
use noodles::core::Region;
pub use noodles::csi;
pub use noodles::tabix;
use noodles::vcf;
use noodles::vcf::record::Chromosome;

pub use noodles_util::variant;

pub mod detect;
use detect::{Compression, Format};
use std::io::Seek;
use std::io::{self, BufRead, BufReader};

pub trait VariantReader {
    fn next_record(&mut self, header: &vcf::Header, v: &mut vcf::Record) -> io::Result<usize>;
}

impl<R> VariantReader for vcf::Reader<R>
where
    R: BufRead,
{
    fn next_record(&mut self, header: &vcf::Header, v: &mut vcf::Record) -> io::Result<usize> {
        self.read_record(header, v)
    }
}
impl<R> VariantReader for bcf::Reader<R>
where
    R: BufRead,
{
    fn next_record(&mut self, header: &vcf::Header, v: &mut vcf::Record) -> io::Result<usize> {
        self.read_record(header, v)
    }
}

pub enum XCF<R> {
    Vcf(Box<dyn VariantReader>),
    IndexedVcf(vcf::IndexedReader<BufReader<Box<R>>>),
    IndexedBcf(bcf::IndexedReader<bgzf::Reader<BufReader<Box<R>>>>),
}

pub struct Reader<R> {
    inner: XCF<R>,
    header: vcf::Header,
    variant: Option<vcf::Record>,
}

impl<R> Reader<R>
where
    R: BufRead + 'static,
{
    pub fn new(inner: XCF<R>, header: vcf::Header) -> Self {
        Self {
            inner,
            header,
            variant: None,
        }
    }

    pub fn next_record(&mut self, header: &vcf::Header, v: &mut vcf::Record) -> io::Result<usize> {
        // if self.variant is set, then use that to set v
        if let Some(variant) = self.variant.take() {
            *v = variant;
            return Ok(1);
        }
        match &mut self.inner {
            XCF::Vcf(reader) => reader.next_record(header, v),
            XCF::IndexedVcf(reader) => reader.read_record(header, v),
            XCF::IndexedBcf(reader) => reader.read_record(header, v),
        }
    }

    pub fn from_reader(reader: Box<R>, path: Option<&str>) -> io::Result<Reader<R>> {
        let mut reader = BufReader::new(reader);
        // this is clearly only available if the reader has type BufRead Seek.
        eprintln!("detecting compression and format");
        let compression = detect::detect_compression(&mut reader)?;
        eprintln!("detecting format");
        let format = detect::detect_format(&mut reader, compression)?;
        let csi = find_index(path);
        eprintln!("csi: {:?}", csi);

        Ok(match (format, compression, csi /*, seekable */) {
            (Format::Vcf, None, _) => {
                let mut reader = vcf::Reader::new(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(Box::new(reader)), header)
            }
            (Format::Vcf, Some(Compression::Bgzf), Some(csi)) => {
                let mut reader = vcf::IndexedReader::new(reader, csi);
                let header = reader.read_header()?;
                Reader::new(XCF::IndexedVcf(reader), header)
            }
            (Format::Vcf, Some(Compression::Bgzf), None) => {
                let bgzf_reader = bgzf::Reader::new(reader);
                let mut reader = vcf::Reader::new(bgzf_reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(Box::new(reader)), header)
            }

            (Format::Bcf, Some(Compression::Bgzf), Some(csi)) => {
                //let bgzf_reader = bgzf::Reader::new(reader);
                let mut reader = bcf::IndexedReader::new(reader, csi);
                let header = reader.read_header()?;
                Reader::new(XCF::IndexedBcf(reader), header)
            }
            (Format::Bcf, _, _) => {
                let mut reader = bcf::Reader::new(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(Box::new(reader)), header)
            }
        })
    }
    pub fn header(&mut self) -> &vcf::Header {
        &self.header
    }
}

impl<R> Reader<R>
where
    R: BufRead + Seek + 'static,
{
    pub fn skip_to(&mut self, header: &vcf::Header, region: &Region) -> io::Result<()> {
        match &mut self.inner {
            XCF::IndexedVcf(reader) => {
                let mut q = reader.query(header, region)?;
                self.variant = match q.next() {
                    Some(Ok(v)) => Some(v),
                    Some(Err(e)) => return Err(e),
                    None => None,
                };
                Ok(())
            }
            _ => unimplemented!(),
        }
    }
}

#[inline]
pub fn chrom_equals(c: &Chromosome, name: &str) -> bool {
    match c {
        Chromosome::Name(n) => n == name,
        Chromosome::Symbol(_) => false,
    }
}

fn find_index(path: Option<&str>) -> Option<csi::Index> {
    if let Some(path) = path {
        let mut index_path = path.clone().to_string();
        // TODO: look for csi or .tbi or use ##idx## in path.
        index_path.push_str(".csi");
        csi::read(index_path)
            .or_else(|_| {
                let mut index_path = path.clone().to_string();
                index_path.push_str(".tbi");
                tabix::read(index_path)
            })
            .ok()
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use noodles::core::{Position, Region};

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
        let _cursor = Cursor::new(vcf_data);
        let path = "tests/t.vcf.gz";
        let rdr = BufReader::new(std::fs::File::open(&path).unwrap());
        let rdr = Box::new(rdr);

        let mut reader = Reader::from_reader(rdr, Some(path)).expect("error creating new reader");
        let header = reader.header().clone();

        let start = Position::try_from(2000).expect("error creating start");
        let stop = Position::try_from(2100).expect("error creating stop");
        let region = Region::new("chr1", start..=stop);

        reader
            .skip_to(&header, &region)
            .expect("error skipping to region");

        let mut v = vcf::Record::default();

        reader
            .next_record(&header, &mut v)
            .expect("error reading record");
        eprintln!("v: {}", v);
    }
}

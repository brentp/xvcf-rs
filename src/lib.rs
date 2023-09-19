pub use noodles::bcf;
use noodles::bgzf;
pub use noodles::csi;
pub use noodles::tabix;
use noodles::vcf;
use noodles::vcf::record::Chromosome;

pub use noodles_util::variant;

pub mod detect;
use detect::{Compression, Format};
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

pub enum XCF {
    Vcf(Box<dyn VariantReader>),
    IndexedVcf(vcf::IndexedReader<BufReader<Box<dyn BufRead>>>),
    IndexedBcf(bcf::IndexedReader<bgzf::Reader<BufReader<Box<dyn BufRead>>>>),
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

    pub fn from_reader(reader: Box<dyn BufRead>, path: Option<&str>) -> io::Result<Reader> {
        let mut reader = BufReader::new(reader);
        eprintln!("detecting compression and format");
        let compression = detect::detect_compression(&mut reader)?;
        eprintln!("detecting format");
        let format = detect::detect_format(&mut reader, compression)?;
        let csi = find_index(path);
        eprintln!("csi: {:?}", csi);

        Ok(match (format, compression, csi) {
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

#[inline]
fn chrom_equals(c: &Chromosome, name: &str) -> bool {
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
        let cursor = Cursor::new(vcf_data);
        let path = "tests/t.vcf";
        let mut rdr = BufReader::new(std::fs::File::open(&path).unwrap());
        let rdr = Box::new(rdr);

        let mut reader = Reader::from_reader(rdr, Some(path)).expect("error creating new reader");
        let header = reader.header().clone();

        let start = Position::try_from(2000).expect("error creating start");
        let stop = Position::try_from(2100).expect("error creating stop");
        let region = Region::new("chr1", start..=stop);

        let mut v = vcf::Record::default();
    }
}

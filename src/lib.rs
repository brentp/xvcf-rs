pub use noodles::bcf;
use noodles::bgzf;
pub use noodles::csi;
use noodles::vcf::record::Chromosome;
pub use noodles::vcf::{self};
use noodles::vcf::{IndexedReader, VariantReader};

pub use noodles_util::variant;

pub mod detect;
use detect::{Compression, Format};
use std::io::{self, Seek};
use std::io::{BufRead, BufReader};

trait BufReadSeek: BufRead + Seek {}

pub enum XCF<R> {
    Vcf(Box<dyn VariantReader<R>>),
    IndexedVcf(vcf::IndexedReader<bgzf::Reader<R>>),
    IndexedBcf(bcf::IndexedReader<bgzf::Reader<R>>),
}

pub struct Reader<R> {
    inner: XCF<R>,
    header: vcf::Header,
    variant: Option<vcf::Record>,
}

impl<R> Reader<R>
where
    R: BufRead,
    vcf::Reader<BufReader<Box<dyn BufRead>>>: VariantReader<R>,
{
    pub fn new(inner: XCF<R>, header: vcf::Header) -> Self {
        Self {
            inner,
            header,
            variant: None,
        }
    }

    pub fn from_reader(reader: Box<dyn BufRead>, path: Option<String>) -> io::Result<Reader<R>> {
        let mut reader = BufReader::new(reader);
        let compression = detect::detect_compression(&mut reader)?;
        let format = detect::detect_format(&mut reader, compression)?;
        let csi = find_index(path);

        let mut rdr: Reader<R> = match (format, compression, csi) {
            (Format::Vcf, None, _) => {
                let mut reader = vcf::Reader::new(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(Box::new(reader)), header)
            }
            (Format::Vcf, Some(Compression::Bgzf), Some(csi)) => {
                let mut bgzf_reader = bgzf::Reader::new(reader);
                let mut reader = IndexedReader::new(bgzf_reader, csi);
                let header = reader.read_header()?;
                Reader::new(XCF::IndexedVcf(reader), header)
            }
            _ => unimplemented!(),
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
            XCF::IndexedVcf(reader) => reader.read_record(header, v),
            XCF::IndexedBcf(reader) => reader.read_record(header, v),
        }
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

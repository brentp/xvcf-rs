use indexmap::IndexMap;
pub use noodles::bcf;
use noodles::bgzf;
use noodles::core::{Position, Region};
pub use noodles::csi;
pub use noodles::tabix;
use noodles::vcf;
use noodles::vcf::header::record::value::map::contig::Name;
use noodles::vcf::header::record::value::map::Contig;
use noodles::vcf::header::record::value::Map;
use noodles::vcf::record::Chromosome;

pub use noodles_util::variant;

pub mod detect;
use detect::{Compression, Format};
use std::io::{self, BufRead, BufReader};
use std::io::{Read, Seek};

pub trait VariantReader: Send {
    fn next_record(&mut self, header: &vcf::Header, v: &mut vcf::Record) -> io::Result<usize>;
}

impl<R> VariantReader for vcf::Reader<R>
where
    R: BufRead + Send,
{
    fn next_record(&mut self, header: &vcf::Header, v: &mut vcf::Record) -> io::Result<usize> {
        self.read_record(header, v)
    }
}
impl<R> VariantReader for bcf::Reader<R>
where
    R: BufRead + Send,
{
    fn next_record(&mut self, header: &vcf::Header, v: &mut vcf::Record) -> io::Result<usize> {
        self.read_record(header, v)
    }
}

pub trait ReadSeek: Read + Seek + 'static + Send {}
pub trait ReadSend: Read + Send {}

pub enum XCF {
    Vcf(Box<dyn VariantReader>),
    IndexedVcf(vcf::IndexedReader<BufReader<Box<dyn ReadSeek>>>),
    IndexedBcf(bcf::IndexedReader<bgzf::Reader<BufReader<Box<dyn ReadSeek>>>>),
}

unsafe impl Send for XCF {}

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

    pub fn from_seek_reader(reader: Box<dyn ReadSeek>, path: Option<&str>) -> io::Result<Reader> {
        let mut breader = BufReader::new(reader);
        let compression = detect::detect_compression(&mut breader)?;
        let format = detect::detect_format(&mut breader, compression)?;
        let csi = find_index(path);
        match (format, compression, csi /*, seekable */) {
            (Format::Bcf, Some(Compression::Bgzf), Some(csi)) => {
                //let bgzf_reader = bgzf::Reader::new(reader);
                let mut reader = bcf::IndexedReader::new(breader, csi);
                let header = reader.read_header()?;
                Ok(Reader::new(XCF::IndexedBcf(reader), header))
            }
            (Format::Vcf, Some(Compression::Bgzf), Some(csi)) => {
                //let bgzf_reader = bgzf::Reader::new(reader);
                let mut reader = vcf::IndexedReader::new(breader, csi);
                let header = reader.read_header()?;
                Ok(Reader::new(XCF::IndexedVcf(reader), header))
            }
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "unsupported format or lack of index",
            )),
        }
    }

    pub fn from_reader<R: Read + Send + 'static>(
        reader: R,
        path: Option<&str>,
    ) -> io::Result<Reader> {
        let mut reader = BufReader::new(reader);
        let compression = detect::detect_compression(&mut reader)?;
        let format = detect::detect_format(&mut reader, compression)?;
        let csi = find_index(path);

        Ok(match (format, compression, csi /*, seekable */) {
            (Format::Vcf, None, _) => {
                let mut reader = vcf::Reader::new(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(Box::new(reader)), header)
            }
            (Format::Vcf, Some(Compression::Bgzf), _) => {
                let bgzf_reader = bgzf::Reader::new(reader);
                let mut reader = vcf::Reader::new(bgzf_reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(Box::new(reader)), header)
            }
            (Format::Bcf, _, _) => {
                let mut reader = bcf::Reader::new(reader);
                let header = reader.read_header()?;
                Reader::new(XCF::Vcf(Box::new(reader)), header)
            }
        })
    }
    pub fn header(&self) -> &vcf::Header {
        &self.header
    }
}

// if the vcf header has the correct info, we check that the contig is present and
// the the region is not beyond the end of the contig.
// if the header does not have this info, then we do not check the contig.
fn check_contig(contigs: &IndexMap<Name, Map<Contig>>, region: &Region) -> io::Result<usize> {
    let contig_i = if !contigs.is_empty() {
        match contigs.get_index_of(region.name()) {
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown contig: {}", region.name()),
                ))
            }
            Some(i) => {
                let contig = contigs.get_index(i).unwrap().1;
                if let Some(len) = contig.length() {
                    if len < usize::from(region.interval().start().unwrap_or(Position::MAX)) {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("region {} is out of bounds for contig {:?}", region, contig),
                        ));
                    }
                }

                i
            }
        }
    } else {
        usize::MAX
    };
    Ok(contig_i)
}

fn simple_skip(reader: &mut Reader, header: &vcf::Header, region: &Region) -> io::Result<()> {
    let mut v = vcf::Record::default();
    let one = Position::try_from(1).unwrap();
    let start =
        vcf::record::Position::try_from(usize::from(region.interval().start().unwrap_or(one)))
            .unwrap();

    let contigs = header.contigs();
    let contig_i = check_contig(contigs, region)?; // check that the contig is valid (if it exists)

    let mut last_chrom: Option<Chromosome> = None;

    loop {
        reader.next_record(header, &mut v)?;
        let end = match v.end() {
            Ok(p) => p,
            _ => vcf::record::Position::from(usize::from(v.position()) + 1),
        };
        // TODO: not sure about 1-based vs 0-based in noodles. might need end > start
        if end >= start && chrom_equals(v.chromosome(), region.name()) {
            reader.variant = Some(v);
            break;
        }

        if contig_i == usize::MAX {
            continue;
        }
        // we check that user isn't requesting a contig that is before the current one.
        if let Some(ref ilast_chrom) = last_chrom {
            // if we have just seen this chromosome, then we can skip the check
            if ilast_chrom == v.chromosome() {
                continue;
            }
            match contigs.get_index_of(match v.chromosome() {
                Chromosome::Symbol(_) => continue,
                Chromosome::Name(n) => n.as_str(),
            }) {
                None => continue,
                Some(i) => {
                    if i > contig_i {
                        // return an error about ordering of chromosomes
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!(
                                "contig {} is out of order relative to {}",
                                v.chromosome(),
                                region.name()
                            ),
                        ));
                    }
                }
            }
            last_chrom = Some(v.chromosome().clone());
        } else {
            last_chrom = Some(v.chromosome().clone());
        }
    }
    Ok(())
}

trait Skip {
    fn skip_to(&mut self, header: &vcf::Header, region: &Region) -> io::Result<()>;
}

impl Skip for Reader {
    fn skip_to(&mut self, header: &vcf::Header, region: &Region) -> io::Result<()> {
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
            XCF::IndexedBcf(reader) => {
                let mut q = reader.query(header, region)?;
                self.variant = match q.next() {
                    Some(Ok(v)) => Some(v),
                    Some(Err(e)) => return Err(e),
                    None => None,
                };
                Ok(())
            }

            _ => simple_skip(self, header, region),
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
        let mut index_path = path.to_string();
        // TODO: look for csi or .tbi or use ##idx## in path.
        index_path.push_str(".csi");
        csi::read(index_path)
            .or_else(|_| {
                let mut index_path = path.to_string();
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
        let rdr = std::fs::File::open(&path).unwrap();
        //let rdr = _cursor;
        let rdr = Box::new(rdr) as Box<dyn Read + Send + 'static>;

        let mut reader = Reader::from_reader(rdr, Some(path)).expect("error creating new reader");
        let header = reader.header().clone();

        let start = Position::try_from(2000).expect("error creating start");
        let region = Region::new("chr1", start..);

        reader
            .skip_to(&header, &region)
            .expect("error skipping to region");

        let mut v = vcf::Record::default();

        reader
            .next_record(&header, &mut v)
            .expect("error reading record");
        assert_eq!(v.position(), Position::try_from(2000).unwrap());
        eprintln!("v: {}", v);
    }
}

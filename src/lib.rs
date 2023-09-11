pub use noodles::bcf;
use noodles::bgzf;
use noodles::core::{Position, Region};
pub use noodles::csi;
pub use noodles::vcf::{self};

pub use noodles_util::variant;

pub mod detect;
use detect::{Compression, Format};
use std::io::{self, Read, Seek};
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

impl<R> Reader<R>
where
    R: BufRead + Seek,
{
    // query simply sets the file pointer to the start of region.
    pub fn query(&mut self, header: &vcf::Header, region: &Region) -> io::Result<()> {
        let p = usize::from(
            region
                .interval()
                .start()
                .unwrap_or(Position::try_from(1).unwrap()),
        );
        loop {
            let mut record = vcf::Record::default();
            self.read_record(header, &mut record)?;
            let var_end = usize::from(
                record
                    .end()
                    .unwrap_or(vcf::record::Position::try_from(usize::MAX).unwrap()),
            );
            if var_end >= p {
                self.variant = Some(record);
                break;
            }
        }

        let tid = header
            .contigs()
            .get(region.name())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid contig"))?;
        let tid = tid.idx().ok_or(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid contig",
        ))?;

        let index = self.index.as_ref().ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "index required for query")
        })?;
        let chunks = index
            .query(tid, region.interval())
            .or_else(|_| Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid query")))?;
        let mut buf = Vec::new();
        let chunk = chunks.iter().next().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid query: no chunks found",
            )
        })?;
        let index = self.index.unwrap();

        match &mut self.inner {
            XCF::Vcf(reader) => reader.seek(chunk.start())?,
            XCF::Bcf(reader) => reader.seek(chunk.start())?,
            XCF::IndexedVcf(reader) => reader.seek(chunk.start())?,
            XCF::IndexedBcf(reader) => reader.seek(chunk.start())?,
        }
        Ok(())
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

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}

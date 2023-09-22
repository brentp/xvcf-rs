/// NOTE!!! these are taken from noodles-util
/// by Michael Macias under MIT license
use std::io::{self, BufRead, Read};

/// A variant format.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    /// Variant Call Format (VCF).
    Vcf,
    /// BCF.
    Bcf,
}

/// A variant compression.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Compression {
    /// BGZF compression.
    Bgzf,
}

pub(crate) fn detect_compression<R>(reader: &mut R) -> io::Result<Option<Compression>>
where
    R: BufRead,
{
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let src = reader.fill_buf()?;

    if let Some(buf) = src.get(..GZIP_MAGIC_NUMBER.len()) {
        if buf == GZIP_MAGIC_NUMBER {
            return Ok(Some(Compression::Bgzf));
        }
    }

    Ok(None)
}

pub(crate) fn detect_format<R>(
    reader: &mut R,
    compression: Option<Compression>,
) -> io::Result<Format>
where
    R: BufRead,
{
    use flate2::bufread::MultiGzDecoder;

    const BCF_MAGIC_NUMBER: [u8; 3] = *b"BCF";
    const VCF_HEADER: [u8; 16] = *b"##fileformat=VCF";

    let src = reader.fill_buf()?;

    if let Some(compression) = compression {
        if compression == Compression::Bgzf {
            let mut decoder = MultiGzDecoder::new(src);
            let mut buf = [0; BCF_MAGIC_NUMBER.len()];
            decoder.read_exact(&mut buf)?;

            if buf == BCF_MAGIC_NUMBER {
                return Ok(Format::Bcf);
            }
            // check that the file is a VCF file. should start with ##fileformat=VCF
            let mut buf = [0; VCF_HEADER.len()];
            decoder.read_exact(&mut buf)?;
            if buf == VCF_HEADER {
                return Ok(Format::Vcf);
            }
            // return error about unknown format
            return Err(io::Error::new(io::ErrorKind::InvalidData, "unknown format"));
        }
        // return error about compression not supported
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unsupported compression",
        ));
    } else if let Some(buf) = src.get(..BCF_MAGIC_NUMBER.len()) {
        if buf == BCF_MAGIC_NUMBER {
            return Ok(Format::Bcf);
        }
    }
    // check for vcf format
    if let Some(buf) = src.get(..VCF_HEADER.len()) {
        if buf == VCF_HEADER {
            return Ok(Format::Vcf);
        }
    }

    // unknown format
    Err(io::Error::new(io::ErrorKind::InvalidData, "unknown format"))
}

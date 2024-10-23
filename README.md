[![Continuous integration](https://github.com/brentp/xvcf-rs/actions/workflows/main.yml/badge.svg)](https://github.com/brentp/xvcf-rs/actions/workflows/main.yml)

This is a small library that attempts to allow efficient reading
of VCF and BCF files that are either compressed or uncompressed and indexed or not.
[rust_htslib](https://github.com/rust-bio/rust-htslib) is used for the parsing, this unifies handling above that.

Even when not indexed, it allows "skipping" via iterating until the requested region is reached.
This is useful for [bedder-rs](https://github.com/quinlan-lab/bedder-rs) but might also be useful elsewhere.

The skipping will be most efficient when the files are compressed and indexed.


## NOTES

This currently will may work as expected when there are multiple variants with the same position.
That will depend on how it is called. This will be fixed in the future.

### Implementation

As implemented, only compressed, indexed VCF or BCF will have the absolute highest performance.
Others will require skipping via iterating over the file until the given region is reached.
This assumes that any indexed file will also support `Seek`.
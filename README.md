This is an attempt to better understand some rust properties.

We'd like to have a reader that supports vcf,vcf.gz,bcf, uncompressed bcf with both
seekable and unseekable readers.

AFAICT, the way to do this is with [specialization](https://std-dev-guide.rust-lang.org/policy/specialization.html)

If the file is not indexed and `Seek`able, then we just iterate over the records to skip to a given region.
If it is `Seek`able and has an index, then we use the `query` functionality.

The document indicates that :
> Only specialization using the min_specialization feature should be used.
> The full specialization feature is known to be unsound.

But I can't get this to compile with `min_specialization`.

Have a look [here](https://github.com/brentp/xvcf-rs/blob/main/src/lib.rs)



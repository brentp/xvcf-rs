[![Continuous integration](https://github.com/brentp/xvcf-rs/actions/workflows/main.yml/badge.svg)](https://github.com/brentp/xvcf-rs/actions/workflows/main.yml)

This is an attempt to better understand some rust properties.

We'd like to have a reader that supports vcf,vcf.gz,bcf, uncompressed bcf with both
seekable and unseekable readers.

AFAICT, the way to do this is with [specialization](https://std-dev-guide.rust-lang.org/policy/specialization.html)

If the file is not indexed and `Seek`able, then we just iterate over the records to skip to a given region.
If it is `Seek`able and has an index, then we use the `query` functionality.

Here is the (working) code as implemented: https://github.com/brentp/xvcf-rs/blob/e027632536b6d5f0da8d7b33578590dd2ce25de1/src/lib.rs#L149-L190

The rust documentation indicates that :
> Only specialization using the min_specialization feature should be used.
> The full specialization feature is known to be unsound.

But I can't get this to compile with `min_specialization`. It also seems that even `min_specialization` will
require `nightly` for the foreseeable future.

Have a look [here](https://github.com/brentp/xvcf-rs/blob/main/src/lib.rs)

Perhaps there's completely different, and better way to do this. Let me know.

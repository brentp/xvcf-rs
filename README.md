[![Continuous integration](https://github.com/brentp/xvcf-rs/actions/workflows/main.yml/badge.svg)](https://github.com/brentp/xvcf-rs/actions/workflows/main.yml)

This is a small library that attempts to allow efficient reading
of VCF and BCF files that are either compressed or uncompressed and indexed or not.
[noodles](https://github.com/zaleus/noodles) is used for the parsing, this unifies handling above that.

Even when not indexed, it allows "skipping" via iterating until the requested region is reached.
This is useful for [bedder-rs](https://github.com/quinlan-lab/bedder-rs) but might also be useful elsewhere.

The skipping will be most efficient when the files are compressed and indexed.


## NOTES

### Implementation

As implemented, only compressed, indexed VCF or BCF will have the absolute highest performance.
Others will require skipping via iterating over the file until the given region is reached.
This assumes that any indexed file will also support `Seek`.

### Specialization (abandoned)

This is not used so we can avoid dependence on nightly rust.

AFAICT, one way to do this is with [specialization](https://std-dev-guide.rust-lang.org/policy/specialization.html)

If the file is not indexed and `Seek`able, then we just iterate over the records to skip to a given region.
If it is `Seek`able and has an index, then we use the `query` functionality.

The rust documentation indicates that :
> Only specialization using the min_specialization feature should be used.
> The full specialization feature is known to be unsound.

But I can't get this to compile with `min_specialization`. It also seems that even `min_specialization` will
require `nightly` for the foreseeable future.

Have a look [here](https://github.com/brentp/xvcf-rs/blob/main/src/lib.rs)

Perhaps there's completely different, and better way to do this. Let me know.

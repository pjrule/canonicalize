# canonicalize

This is a tool for exploring the configuration space (or *metagraph*) of districtings of an NxM grid into M districts of equal size via the [recombination move](https://arxiv.org/abs/1911.05725).

## Running
`canonicalize` takes the dimensions of the grid as arguments and reads districting plans from `stdin`. The districting plans should be in the comma-delimited format used by [Zachary Schutzman's `enumerator`](https://github.com/zschutzman/enumerator). Example:

```
cargo run -- --width 7 --height 7 --n-districts 7 < "enum_7,7_7_7_rc.txt"
```

where `enum_7,7_7_7_rc.txt` contains an enumeration of the 7x7 → 7 districting plans.

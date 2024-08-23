# multilinear-prover

Multilinear prover project

## NOTES

Currently, code is written slightly naively and isn't written with an eye for heavy optimization. However, it seems already with just the compiler itself and some tasteful application of simple `rayon` methods, it's getting quite fast. Easy ways to cut corners are:

- A more efficient (especially cache-efficient) FFT
- Costs incurred by clones and copies which happen during the IOP and any sort of extension field lifting

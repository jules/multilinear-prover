# multilinear-prover

Multilinear prover project

## NOTES

Currently, code is written slightly naively and isn't written with an eye for heavy optimization. However, it seems already with just the compiler itself and some tasteful application of simple `rayon` methods, it's getting quite fast. Easy ways to cut corners are:

- A more efficient (especially cache-efficient) FFT
- Costs incurred by clones and copies which happen during the IOP and any sort of extension field lifting

## LEFT TODO

- [x] Prodcheck argument
- [ ] Some lookup argument (likely best choice here is log-derivative)
- [x] Virtual poly toolkit that does more than just N degree mult
- [ ] Basefold (half done implementation somewhere in boojum that i made and can probably port over)
- [ ] FRI (just since i want to pit a few commitment schemes against each other in recursive settings - note that we may instead do this as a RISC-V program)
- [ ] RISC-V circuit (likely doing one read one write per row, may have to exclude memory ops with this configuration for now and figure that out later)
- [ ] FFT optimization
- [ ] Optimizations to zerocheck IOP (papers already in the comments)

And of course general cleanup and documentation

Also thinking about reusing columns between arguments and saving on commitment costs. but it means we cant batch commit certain columns which might incur more cost - though with a linear commitment scheme it should be fine. bench to be sure

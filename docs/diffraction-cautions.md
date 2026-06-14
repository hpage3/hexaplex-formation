# Diffraction Interpretation Cautions

## d-spacing is not a direct PyMOL distance

A diffraction feature at d ≈ 4.5 Å corresponds to reciprocal-space structure near:

q = 2π / d

It does not imply that the dominant contributing atom pairs are separated by a literal 4.5 Å real-space distance.

## Intensities are not additive by atom group

For the full Hexaplex:

F(full) = F(hexads) + F(scaffold)

but

I(full) = |F(full)|²

Therefore:

I(full) = I(hexads) + I(scaffold) + cross/interference terms

The hexads-only and scaffold-only models are useful comparative controls, but they should not be interpreted as clean additive decompositions of the full intensity.

## Current conservative interpretation

The d ≈ 4.5 Å feature appears to be a useful signature of the organized Hexaplex scaffold architecture.

The feature differentiates the Hexaplex scaffold from isolated beta-strand / alanine-pair controls, suggesting that the folded/twisted multi-strand scaffold geometry is important.

This does not yet prove a unique formation pathway or stabilizing mechanism.

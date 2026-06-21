# Nick Rotation Range Ambiguity Note

Nick's shared notebook/code includes:

```python
rotations = range(0, 5, 180)
```

Python interprets this as `range(start=0, stop=5, step=180)`, which produces only the single value `0`. It does not mean rotation sampling from 0 degrees to 180 degrees in 5-degree steps.

Based on Nick's accompanying text, the likely intended interpretation is 5-degree rotation sampling across the 0-180 degree interval. For future corrected reproduction of Nick's stated workflow, prefer:

```python
rotations = range(0, 181, 5)
```

This explicitly includes 180 degrees. Because 180 degrees may be symmetry-equivalent to 0 degrees after image spinning/symmetrization, the corrected analysis should either justify including the endpoint or run a small sensitivity check comparing:

```python
range(0, 180, 5)
range(0, 181, 5)
```

This range ambiguity is separate from Asem's earlier non-accumulating rotation correction. The preserved source archive should remain unchanged; corrected reanalysis should document any reinterpretation in analysis scripts and reports.

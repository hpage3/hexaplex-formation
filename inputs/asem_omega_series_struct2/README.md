# Asem Omega-Series Structures from struct2

Imported on 2026-06-21 from:

`C:\Users\hpage3\OneDrive - Georgia Institute of Technology\Documents\GitHub\research\struct2`

These PDB files are Asem's manually constructed omega-series structures derived from `1_558214.pdb`. Asem reported that he progressively increased the peptide dihedral while holding all other atoms fixed. The tradeoff is that this also increases the N-CA(next)-C angle, introducing angle strain that may be worse energetically than the original peptide-dihedral strain. Asem suggested a value around 172 degrees may be acceptable.

This is a controlled one-parameter omega series. The goal is not omega = 180 degrees at all costs; the goal is a best compromise, likely near 172 degrees if H-bonding, sterics, and later diffraction-ready periodic modeling remain acceptable.

| Target omega (deg) | Asem peptide dihedral (deg) | Asem N-CA(next)-C angle (deg) | File |
| ---: | ---: | ---: | --- |
| 167 | 167.00 | 113.01 | `1_558214_omega167.pdb` |
| 168 | 168.00 | 113.81 | `1_558214_omega168.pdb` |
| 169 | 169.00 | 114.62 | `1_558214_omega169.pdb` |
| 170 | 170.00 | 115.43 | `1_558214_omega170.pdb` |
| 171 | 171.00 | 116.25 | `1_558214_omega171.pdb` |
| 172 | 172.00 | 117.08 | `1_558214_omega172.pdb` |
| 173 | 173.00 | 117.91 | `1_558214_omega173.pdb` |
| 174 | 174.00 | 118.75 | `1_558214_omega174.pdb` |
| 175 | 175.00 | 119.59 | `1_558214_omega175.pdb` |
| 176 | 176.00 | 120.43 | `1_558214_omega176.pdb` |
| 177 | 177.00 | 121.28 | `1_558214_omega177.pdb` |
| 178 | 178.00 | 122.13 | `1_558214_omega178.pdb` |
| 179 | 179.00 | 122.99 | `1_558214_omega179.pdb` |
| 180 | 180.00 | 123.85 | `1_558214_omega180.pdb` |

The un-suffixed `1_558214.pdb` is included as the original baseline copy from the same source folder.

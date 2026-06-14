"""Small PDB parsing helpers for structure-level Hexaplex comparisons."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PDBAtom:
    record_type: str
    atom_serial: int | None
    atom_name: str
    alt_loc: str
    residue_name: str
    chain_id: str
    residue_number: int | None
    insertion_code: str
    x: float
    y: float
    z: float
    occupancy: float | None
    temp_factor: float | None
    element: str


def _slice(line: str, start: int, end: int) -> str:
    return line[start:end] if len(line) > start else ""


def _parse_int(value: str) -> int | None:
    value = value.strip()
    return int(value) if value else None


def _parse_float(value: str) -> float | None:
    value = value.strip()
    return float(value) if value else None


def _infer_element(atom_name: str) -> str:
    stripped = atom_name.strip()
    if not stripped:
        return ""

    letters = [char for char in stripped if char.isalpha()]
    if not letters:
        return ""

    # PDB atom names often encode right-justified protein atoms such as " CA "
    # for carbon alpha. Prefer the first alphabetic character for conservative
    # hydrogen filtering and blank-element robustness.
    return letters[0].upper()


def parse_pdb_atom_line(line: str, line_number: int | None = None) -> PDBAtom:
    record_type = _slice(line, 0, 6).strip()
    if record_type not in {"ATOM", "HETATM"}:
        raise ValueError(f"Line is not an ATOM/HETATM record: {line.rstrip()}")

    try:
        x = float(_slice(line, 30, 38))
        y = float(_slice(line, 38, 46))
        z = float(_slice(line, 46, 54))
    except ValueError as exc:
        where = f" on line {line_number}" if line_number is not None else ""
        raise ValueError(f"Invalid PDB coordinates{where}: {line.rstrip()}") from exc

    atom_name = _slice(line, 12, 16).strip()
    element = _slice(line, 76, 78).strip()
    if not element:
        element = _infer_element(atom_name)

    return PDBAtom(
        record_type=record_type,
        atom_serial=_parse_int(_slice(line, 6, 11)),
        atom_name=atom_name,
        alt_loc=_slice(line, 16, 17).strip(),
        residue_name=_slice(line, 17, 20).strip(),
        chain_id=_slice(line, 21, 22).strip(),
        residue_number=_parse_int(_slice(line, 22, 26)),
        insertion_code=_slice(line, 26, 27).strip(),
        x=x,
        y=y,
        z=z,
        occupancy=_parse_float(_slice(line, 54, 60)),
        temp_factor=_parse_float(_slice(line, 60, 66)),
        element=element,
    )


def load_pdb_atoms(path: str | Path) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(("ATOM", "HETATM")):
                atoms.append(parse_pdb_atom_line(line, line_number=line_number))
    return atoms


def is_hydrogen(atom: PDBAtom) -> bool:
    return atom.element.upper() == "H" or atom.atom_name.strip().upper().startswith("H")


def heavy_atoms(atoms: Iterable[PDBAtom]) -> list[PDBAtom]:
    return [atom for atom in atoms if not is_hydrogen(atom)]


def atom_identity_key(
    atom: PDBAtom,
) -> tuple[str, str, str, str, int | None, str, float, float, float, str]:
    return (
        atom.record_type,
        atom.atom_name,
        atom.residue_name,
        atom.chain_id,
        atom.residue_number,
        atom.insertion_code,
        round(atom.x, 3),
        round(atom.y, 3),
        round(atom.z, 3),
        atom.element.upper(),
    )


def dedupe_exact_atoms(atoms: Iterable[PDBAtom]) -> list[PDBAtom]:
    seen: set[tuple[str, str, str, str, int | None, str, float, float, float, str]] = set()
    deduped: list[PDBAtom] = []
    for atom in atoms:
        key = atom_identity_key(atom)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(atom)
    return deduped


def _format_atom_name(atom_name: str, element: str) -> str:
    atom_name = atom_name[:4]
    if len(atom_name) == 4:
        return atom_name
    if len(element.strip()) == 1:
        return f" {atom_name:<3}"[:4]
    return f"{atom_name:<4}"[:4]


def write_pdb_atoms(atoms: Iterable[PDBAtom], path: str | Path) -> None:
    """Write atoms as simple fixed-width PDB records.

    Output atom serials are renumbered sequentially to keep cleaned structures
    compact and unambiguous after hydrogen filtering and deduplication.
    """

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        for serial, atom in enumerate(atoms, start=1):
            record_type = atom.record_type if atom.record_type in {"ATOM", "HETATM"} else "ATOM"
            residue_number = atom.residue_number if atom.residue_number is not None else 0
            occupancy = atom.occupancy if atom.occupancy is not None else 1.0
            temp_factor = atom.temp_factor if atom.temp_factor is not None else 0.0
            element = (atom.element or _infer_element(atom.atom_name)).strip().upper()[:2]
            line = (
                f"{record_type:<6}{serial:5d} "
                f"{_format_atom_name(atom.atom_name, element)}"
                f"{atom.alt_loc[:1]:1}"
                f"{atom.residue_name[:3]:>3} "
                f"{atom.chain_id[:1]:1}"
                f"{residue_number:4d}"
                f"{atom.insertion_code[:1]:1}   "
                f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
                f"{occupancy:6.2f}{temp_factor:6.2f}"
                f"          {element:>2}\n"
            )
            handle.write(line)
        handle.write("END\n")


def atom_count(atoms: Iterable[PDBAtom]) -> int:
    return sum(1 for _ in atoms)


def residue_count(atoms: Iterable[PDBAtom]) -> int:
    residues = {
        (
            atom.chain_id,
            atom.residue_number,
            atom.insertion_code,
            atom.residue_name,
        )
        for atom in atoms
    }
    return len(residues)


def chain_ids(atoms: Iterable[PDBAtom]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for atom in atoms:
        if atom.chain_id not in seen:
            seen.add(atom.chain_id)
            ordered.append(atom.chain_id)
    return ordered


def residue_names(atoms: Iterable[PDBAtom]) -> list[str]:
    return sorted({atom.residue_name for atom in atoms if atom.residue_name})


def bounding_box(atoms: Iterable[PDBAtom]) -> tuple[float, float, float, float, float, float] | None:
    atom_list = list(atoms)
    if not atom_list:
        return None

    xs = [atom.x for atom in atom_list]
    ys = [atom.y for atom in atom_list]
    zs = [atom.z for atom in atom_list]
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def centroid(atoms: Iterable[PDBAtom]) -> tuple[float, float, float] | None:
    atom_list = list(atoms)
    if not atom_list:
        return None

    count = len(atom_list)
    return (
        sum(atom.x for atom in atom_list) / count,
        sum(atom.y for atom in atom_list) / count,
        sum(atom.z for atom in atom_list) / count,
    )

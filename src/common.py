"""Read AMBER topologies and handle coordinates in XYZ or AMBER restart format.
The scripts also use these helpers to save results and check the inputs copied into each run.
"""

import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr"
).split()
Z = {element: number for number, element in enumerate(ELEMENTS, start=1)}
CHARGE_SCALE = 18.2223


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def text(path, value):
    Path(path).write_text(value, encoding="utf-8", newline="\n")


def write(path, value):
    path = Path(path)
    payload = json.dumps(value, indent=2, allow_nan=False) + "\n"
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=path.parent,
        prefix="." + path.name,
        delete=False,
    ) as f:
        temporary = f.name
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(temporary, path)


def seal(root, names, manifest="manifest.json"):
    root = Path(root)
    write(root / manifest, {str(Path(p).as_posix()): sha(root / p) for p in sorted(names)})


def verify(root, manifest="manifest.json"):
    root = Path(root).resolve()
    for relative, digest in read(root / manifest).items():
        p = (root / relative).resolve()
        if not p.is_relative_to(root) or sha(p) != digest:
            raise ValueError("Frozen file changed: " + relative)


def xyz(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Missing XYZ header in {path}")
    atoms = []
    for line in lines[2:]:
        fields = line.split()
        if fields:
            element, *coordinates = fields
            atoms.append((element, [float(value) for value in coordinates]))
    if len(atoms) != int(lines[0]) or any(len(coords) != 3 or element not in Z for element, coords in atoms):
        raise ValueError("Malformed XYZ: " + str(path))
    if any(not math.isfinite(value) for _, coords in atoms for value in coords):
        raise ValueError("Nonfinite XYZ")
    return atoms


def write_xyz(path, atoms, title):
    lines = [str(len(atoms)), title]
    for element, coordinates in atoms:
        position = " ".join(f"{value:.12f}" for value in coordinates)
        lines.append(f"{element} {position}")
    text(path, "\n".join(lines) + "\n")


def mask(indices, prefix="@"):
    """Convert zero-based indices to an AMBER mask such as @1-3,7."""
    ranges = []
    for index in sorted(set(indices)):
        amber_index = index + 1
        if ranges and amber_index == ranges[-1][1] + 1:
            ranges[-1][1] = amber_index
        else:
            ranges.append([amber_index, amber_index])
    return prefix + ",".join(str(start) if start == end else f"{start}-{end}" for start, end in ranges)


def write_rst(path, coords, title="Reconstructed fixed-QM MM environment"):
    values = [value for position in coords for value in position]
    lines = [title, f"{len(coords):6d}"]
    for start in range(0, len(values), 6):
        lines.append("".join(f"{value:12.7f}" for value in values[start : start + 6]))
    text(path, "\n".join(lines) + "\n")


def read_rst(path, natom):
    lines = Path(path).read_text().splitlines()
    if int(lines[1].split()[0]) != natom:
        raise ValueError("Restart atom count changed")
    values = []
    for line in lines[2:]:
        values.extend(
            float(line[i : i + 12]) for i in range(0, len(line), 12) if line[i : i + 12].strip()
        )
        if len(values) >= natom * 3:
            break
    if len(values) < natom * 3 or not all(map(math.isfinite, values[: natom * 3])):
        raise ValueError("Truncated/nonfinite restart coordinates")
    return [values[i : i + 3] for i in range(0, natom * 3, 3)]


class Prmtop:
    def __init__(self, path):
        self.lines = Path(path).read_text().splitlines(keepends=True)
        self.sections = {}
        starts = [i for i, s in enumerate(self.lines) if s.startswith("%FLAG ")]
        for k, i in enumerate(starts):
            end = starts[k + 1] if k + 1 < len(starts) else len(self.lines)
            fields = self.lines[i].split()
            if len(fields) != 2:
                raise ValueError(f"Invalid prmtop %FLAG on line {i + 1}")
            name = fields[1]
            if name in self.sections:
                raise ValueError("Duplicate prmtop %FLAG " + name)
            j = i + 1
            while j < end and not self.lines[j].startswith("%FORMAT"):
                j += 1
            if j == end:
                raise ValueError("Missing %FORMAT for prmtop %FLAG " + name)
            m = re.fullmatch(
                r"%FORMAT\(\s*(\d*)([aAiIeEfFdD])(\d+)(?:\.(\d+))?\s*\)\s*", self.lines[j]
            )
            if not m or int(m[1] or 1) < 1 or int(m[3]) < 1:
                raise ValueError("Invalid or unsupported %FORMAT for prmtop %FLAG " + name)
            repeat, kind, width, precision = m.groups()
            width = int(width)
            records = [
                s.rstrip("\r\n") for s in self.lines[j + 1 : end] if not s.startswith("%COMMENT")
            ]
            tokens = [s[x : x + width] for s in records for x in range(0, len(s), width)]
            if kind.lower() == "a":
                values = [s.strip() for s in tokens]
                while values and not values[-1]:
                    values.pop()
            elif kind.lower() == "i":
                values = [int(s) for s in tokens if s.strip()]
            else:
                values = [float(s.upper().replace("D", "E")) for s in tokens if s.strip()]
            self.sections[name] = dict(
                start=i,
                end=end,
                values=values,
                kind=kind,
                width=width,
                repeat=int(repeat or 1),
                precision=int(precision or 0),
            )

    def get(self, name):
        return self.sections[name]["values"]

    def bonds(self):
        result = set()
        for section in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
            # Each triple stores two coordinate offsets (3 * atom index) and a parameter index.
            bonds = self.get(section)
            for start in range(0, len(bonds), 3):
                first = bonds[start] // 3
                second = bonds[start + 1] // 3
                result.add(tuple(sorted((first, second))))
        return result


def run_command(argv, cwd, log, env=None):
    with Path(log).open("xb") as f:
        result = subprocess.run(
            list(map(str, argv)), cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT
        )
    if result.returncode:
        raise RuntimeError(f"Command returned {result.returncode}: {argv[0]}; inspect {log}")

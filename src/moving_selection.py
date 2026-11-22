"""Turn a list of selected atoms into an AMBER mask of at most 256 characters.
Build it from the topology's atom names and indices, then check that it selects exactly the requested atoms.
"""

from collections import defaultdict
import re

from common import mask as numeric_mask


LIMIT = 256
_NAME = re.compile(r"[A-Za-z0-9_+]+\Z")


def _topology(parm):
    count = parm.get("POINTERS")[0]
    if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
        raise ValueError("Moving-mask NATOM must be a positive integer")
    names = parm.get("ATOM_NAME")
    starts = parm.get("RESIDUE_POINTER")
    if len(names) != count or not starts or starts[0] != 1:
        raise ValueError("Incomplete atom-name or residue-pointer fields")
    if any(isinstance(i, bool) or not isinstance(i, int) or i < 1 or i > count for i in starts):
        raise ValueError("Invalid residue pointer")
    if starts != sorted(set(starts)):
        raise ValueError("Residue pointers are not strictly increasing")
    residues = []
    residue_for_atom = [None] * count
    for ri, start in enumerate(starts):
        end = starts[ri + 1] - 1 if ri + 1 < len(starts) else count
        indices = list(range(start - 1, end))
        if any(not isinstance(names[i], str) or not names[i] for i in indices):
            raise ValueError("Empty/invalid atom name")
        # Duplicate names cannot describe an arbitrary subset safely. The
        # reconstructed model maps unique canonical names in every residue.
        if len({names[i] for i in indices}) != len(indices):
            raise ValueError("Nonunique atom names within residue " + str(ri + 1))
        residues.append(indices)
        for i in indices:
            residue_for_atom[i] = ri
    return count, names, residues, residue_for_atom


def _numeric(text, limit, label):
    chosen = set()
    if not text:
        raise ValueError("Empty " + label + " selector")
    for part in text.split(","):
        match = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not match:
            raise ValueError("Invalid numerical " + label + " selector: " + part)
        a = int(match[1])
        b = int(match[2] or a)
        if not 1 <= a <= b <= limit:
            raise ValueError("Out-of-range " + label + " selector: " + part)
        chosen.update(range(a - 1, b))
    return chosen


def evaluate_mask(parm, expression):
    """Independently evaluate the restricted subset used by this generator.

    Compound :residue@atom expressions mean residue AND atom selection, as
    documented by AMBER. ! has higher precedence than & and then |.
    This function is not a general-purpose AMBER mask implementation.
    """
    count, names, residues, residue_for_atom = _topology(parm)
    universe = set(range(count))
    tokens = []
    offset = 0
    while offset < len(expression):
        char = expression[offset]
        if char.isspace():
            offset += 1
            continue
        if char in "()&|!":
            tokens.append(char)
            offset += 1
            continue
        if char not in ":@":
            raise ValueError("Unsupported moving-mask character: " + char)
        end = offset + 1
        while end < len(expression) and expression[end] not in ":@()&|! \t\r\n":
            end += 1
        if end == offset + 1:
            raise ValueError("Empty moving-mask selector")
        tokens.append(expression[offset:end])
        offset = end
    position = 0

    def select(token):
        body = token[1:]
        if token[0] == ":":
            chosen = _numeric(body, len(residues), "residue")
            return {i for i, ri in enumerate(residue_for_atom) if ri in chosen}
        if re.fullmatch(r"\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*", body):
            return _numeric(body, count, "atom")
        literals = body.split(",")
        if any(not _NAME.fullmatch(name) for name in literals):
            raise ValueError("Unsupported atom-name selector: " + body)
        wanted = set(literals)
        return {i for i, name in enumerate(names) if name in wanted}

    def primary():
        nonlocal position
        if position >= len(tokens):
            raise ValueError("Unexpected end of moving mask")
        token = tokens[position]
        position += 1
        if token == "(":
            result = expression_or()
            if position >= len(tokens) or tokens[position] != ")":
                raise ValueError("Unbalanced moving-mask parentheses")
            position += 1
            return result
        if token.startswith((":", "@")):
            return select(token)
        raise ValueError("Unexpected moving-mask token: " + token)

    def unary():
        nonlocal position
        if position < len(tokens) and tokens[position] == "!":
            position += 1
            return universe - unary()
        return primary()

    def expression_and():
        nonlocal position
        result = unary()
        while position < len(tokens):
            token = tokens[position]
            if token == "&":
                position += 1
                result &= unary()
            elif token.startswith(("@", ":")):
                # AMBER's compound elementary residue/atom selection.
                result &= unary()
            else:
                break
        return result

    def expression_or():
        nonlocal position
        result = expression_and()
        while position < len(tokens) and tokens[position] == "|":
            position += 1
            result |= expression_and()
        return result

    if not tokens:
        raise ValueError("Empty moving mask")
    selected = expression_or()
    if position != len(tokens):
        raise ValueError("Unexpected trailing moving-mask tokens")
    return sorted(selected)


def build_moving_selection(parm, moving_indices):
    count, names, residues, _ = _topology(parm)
    indices = list(moving_indices)
    if not indices:
        raise ValueError("The movable atom set must not be empty")
    if any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < count for i in indices):
        raise ValueError("Movable atom index is invalid")
    if len(set(indices)) != len(indices):
        raise ValueError("Duplicate movable atom indices")
    target = set(indices)
    candidates = []
    candidates.append(("numeric_atom_indices", numeric_mask(indices, "@")))
    groups = defaultdict(list)
    named_supported = True
    for ri, atoms in enumerate(residues):
        selected = [i for i in atoms if i in target]
        if not selected:
            continue
        if len(selected) == len(atoms):
            key = ("all", ())
        else:
            included = tuple(sorted(names[i] for i in selected))
            excluded = tuple(sorted(names[i] for i in atoms if i not in target))
            if any(not _NAME.fullmatch(name) or name.isdigit() for name in included + excluded):
                named_supported = False
                break
            # Include adds '@'; exclude adds '(&!@)' around the residue
            # selector, four additional characters relative to include.
            if len(",".join(included)) <= len(",".join(excluded)) + 4:
                key = ("include", included)
            else:
                key = ("exclude", excluded)
        groups[key].append(ri)
    if named_supported:
        parts = []
        for (kind, atom_names), residue_indices in groups.items():
            selection = numeric_mask(residue_indices, ":")
            if kind == "include":
                selection += "@" + ",".join(atom_names)
            elif kind == "exclude":
                selection = "(" + selection + "&!@" + ",".join(atom_names) + ")"
            parts.append(selection)
        candidates.append(("grouped_residue_literal_names", "|".join(parts)))
    evidence = []
    for strategy, expression in candidates:
        selected = evaluate_mask(parm, expression)
        if selected != sorted(target):
            raise ValueError("Generated moving mask does not select the exact requested atom set")
        evidence.append(
            dict(
                strategy=strategy,
                mask=expression,
                length=len(expression),
                selected_count=len(selected),
                exact_selection_verified=True,
            )
        )
    viable = [item for item in evidence if item["length"] <= LIMIT]
    if not viable:
        raise ValueError(
            "No exact moving mask fits AMBER 256-character limit; candidate lengths: "
            + str({x["strategy"]: x["length"] for x in evidence})
        )
    chosen = min(viable, key=lambda item: (item["length"], item["mask"]))
    return dict(
        mask=chosen["mask"],
        selected_indices_0based=sorted(target),
        count=len(target),
        length=chosen["length"],
        strategy=chosen["strategy"],
        candidate_masks=evidence,
        exact_selection_verified=True,
        topology_modified=False,
        atom_names_modified=False,
        mask_limit=LIMIT,
        scope="Static exact selection; no distance-dependent membership",
    )

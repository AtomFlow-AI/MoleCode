"""MoleCode stereo extension (profile stereo/1).

Design summary
--------------
MoleCode's base grammar (atoms, bonds, explicit H counts, charge, E/Z) is a
*constitution* language.  Stereo facts do not fit a single fixed granularity:

* a tetrahedral centre is unary and local, but its *absolute* CIP descriptor is a
  global derived quantity (edit far away -> recompute every label);
* enhanced stereo (AND / OR / ABS) and axial chirality are *relational* facts
  that cannot be attached to one atom or one bond;
* isotopes are unary, but they must be applied *before* CIP labels are trusted.

The extension therefore uses four layers:

S0 base projection
    atoms + bonds + H counts + charge + E/Z.  Unchanged, readable by the stock
    MoleCode parser.  It must stay *sound*: it may never assert a stereo fact
    that the extension layer does not entail (in particular, AND/OR members are
    written with **no** absolute _R/_S suffix).

S1 local configuration records (source of truth)
    Pfx_TET_1{tet=Pfx_C_2; par=cw} -- tetrahedral parity of one centre.

S2 relational / higher-order objects
    Pfx_SG_1{sg=AND:1; mem=Pfx_C_2,Pfx_C_4}  enhanced stereo group
    Pfx_AX_1{ax=cw; axis=Pfx_C_7,Pfx_C_8; ref1=Pfx_C_2; ref2=Pfx_C_9} axial
    Pfx_ISO_1{iso=13; atom=Pfx_C_1} isotope (unary, but needs a slot)

S3 derived descriptor cache
    the _R/_S/_r/_s suffix on an atom id (modern CIP, declared in the header as
    cip-impl) and the ===|E| / ===|Z| bond annotation.  These are *claims* that
    the reader re-derives and checks; they are never the source of truth when S1
    is present.

Reference frame for par
-----------------------
par=cw / par=ccw is the RDKit tetrahedral chiral tag of the atom with respect to
the **document order of its incident edges** (the usual RDKit / SMILES
convention: first neighbour in bond order, implicit H last).  Because MoleCode
writes one bond per line and the reader adds bonds in document order, this frame
is explicit, reproducible and round-trip stable.  New bonds for an existing
centre should therefore be appended (documented), and the CIP cache lets a
validator detect a broken frame.

Invariants (checked by validate_document)
-----------------------------------------
I1 soundness  - no absolute suffix on AND/OR members in the base projection.
I2 descriptor - every suffix is re-derived with the declared CIP implementation.
I3 groups     - each stereo element belongs to at most one group; ABS members
                are not simultaneously AND/OR members.
I4 isotopy    - isotopes are applied before CIP validation.
I5 axial refs - ref atoms must be neighbours of the two axis atoms.
I6 judge      - full InChIKey is not a valid judge when AND/OR or axial
                chirality is present (InChI drops / is blind to both).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from rdkit import Chem

try:  # modern (Hanson 2018) CIP labeler, available since RDKit 2021.09
    from rdkit.Chem import rdCIPLabeler as _rdCIPLabeler
except Exception:  # pragma: no cover - very old RDKit
    _rdCIPLabeler = None

PROFILE_ID = "stereo/1"
CIP_IMPL = "hanson2018"
CIP_IMPL_DETAIL = "rdkit.Chem.rdCIPLabeler"
PROFILE_HEADER = (
    f"%% molecode-profile: {PROFILE_ID}",
    f"%% cip-impl: {CIP_IMPL} ({CIP_IMPL_DETAIL})",
)
STEREO_SECTION_COMMENT = "%% ---- stereo layer (stereo/1) ----"

PARITY_MODES = ("all", "group_only", "off")

_STEREO_GROUP_TYPE = getattr(Chem, "StereoGroupType", None)
GROUP_TYPE_TO_NAME = {}
if _STEREO_GROUP_TYPE is not None:  # pragma: no branch - always true on supported RDKit
    for _name, _attr in (
        ("AND", "STEREO_AND"),
        ("OR", "STEREO_OR"),
        ("ABS", "STEREO_ABSOLUTE"),
    ):
        _value = getattr(_STEREO_GROUP_TYPE, _attr, None)
        if _value is not None:
            GROUP_TYPE_TO_NAME[_value] = _name
GROUP_NAME_TO_TYPE = {v: k for k, v in GROUP_TYPE_TO_NAME.items()}

_CW = Chem.ChiralType.CHI_TETRAHEDRAL_CW
_CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CCW


def _atrop_cw_ccw():
    """Return the RDKit atropisomer bond-stereo enums (None on old RDKit)."""
    return (
        getattr(Chem.BondStereo, "STEREOATROPCW", None),
        getattr(Chem.BondStereo, "STEREOATROPCCW", None),
    )


def atrop_supported() -> bool:
    """Whether this RDKit build can represent axial (atropisomer) chirality."""
    cw, ccw = _atrop_cw_ccw()
    return cw is not None and ccw is not None


_ATROP = tuple(v for v in _atrop_cw_ccw() if v is not None)

# Atom-id suffixes that the extension understands (R/S = absolute, r/s = pseudo-
# asymmetric).  The stock parser only accepts _R/_S; _r/_s is invisible to it
# (the atom keeps a longer id but the constitution is unchanged -> sound).
_SUFFIX_RE = re.compile(r"^(?P<base>.+)_(?P<cip>[RSrs])$")
_OBJECT_RE = re.compile(r"^\s*(?P<oid>[A-Za-z_][\w]*)\{(?P<body>[^{}]*)\}\s*$")
_MEMBER_EDGE_RE = re.compile(r"^\s*(?P<src>[\w]+)\s*-\.->\s*(?P<dst>[\w]+)\s*$")


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def tag_to_par(tag) -> Optional[str]:
    """RDKit chiral tag -> cw / ccw (document-edge reference frame)."""
    if tag == _CW:
        return "cw"
    if tag == _CCW:
        return "ccw"
    return None


def par_to_tag(par: str):
    p = par.strip().lower()
    if p == "cw":
        return _CW
    if p == "ccw":
        return _CCW
    raise ValueError(f"invalid tetrahedral parity {par!r} (expected cw/ccw)")


def is_lowercase_cip(cip: str) -> bool:
    return cip in ("r", "s")


def group_type_name(gt) -> str:
    return GROUP_TYPE_TO_NAME.get(gt, str(gt))


def prepare_mol(mol: Chem.Mol) -> Chem.Mol:
    """Copy mol and run legacy stereo perception (AssignStereochemistry).

    The modern CIP labeler is deliberately NOT run on the working mol: RDKit's
    rdCIPLabeler rewrites E/Z bond stereo into CIS/TRANS, which would destroy
    the E/Z information needed by the base projection.  Modern CIP descriptors
    are computed on a throw-away copy by modern_cip_map().
    """
    if mol is None:
        raise ValueError("prepare_mol() got None")
    m = Chem.Mol(mol)
    try:
        # cleanIt=False: a cleanIt pass can drop E/Z labels when bond
        # directions disagree (fused dienes share a single bond direction),
        # which would make the base projection lossy on re-serialization.
        Chem.AssignStereochemistry(m, cleanIt=False, force=True)
    except Exception:
        pass
    return m


def modern_cip_map(mol: Chem.Mol) -> Dict[int, str]:
    """Modern (Hanson 2018) CIP labels, computed on a throw-away copy.

    The RDKit labeler also rewrites double bond E/Z into CIS/TRANS as a side
    effect, so it must never run on the mol that is subsequently serialized
    with the base E/Z grammar.
    """
    if mol is None:
        return {}
    m = Chem.Mol(mol)
    try:
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    except Exception:
        pass
    if _rdCIPLabeler is not None:
        try:
            _rdCIPLabeler.AssignCIPLabels(m)
        except Exception:
            pass
    return cip_map(m)


def cip_map(mol: Chem.Mol) -> Dict[int, str]:
    return {
        a.GetIdx(): a.GetProp("_CIPCode")
        for a in mol.GetAtoms()
        if a.HasProp("_CIPCode")
    }


def legacy_cip_map(mol: Chem.Mol) -> Dict[int, str]:
    m = Chem.Mol(mol)
    try:
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    except Exception:
        return {}
    return cip_map(m)


def inchi_key(mol: Chem.Mol) -> str:
    try:
        inchi = Chem.MolToInchi(mol)
        return Chem.InchiToInchiKey(inchi) if inchi else ""
    except Exception:
        return ""


def canonical_smiles(mol: Chem.Mol) -> str:
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return ""


def cx_smiles(mol: Chem.Mol) -> str:
    try:
        return Chem.MolToCXSmiles(mol)
    except Exception:
        return ""


def _neighbors_excluding(mol: Chem.Mol, atom_idx: int, exclude: int) -> List[int]:
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
        if n.GetIdx() != exclude
    ]


def _cx_axial_tokens(mol: Chem.Mol) -> set:
    """Return the set of wU: / wD: CXSMILES tokens of mol."""
    try:
        text = Chem.MolToCXSmiles(mol)
    except Exception:
        return set()
    if "|" not in text:
        return set()
    start = text.index("|")
    ext = text[start:].strip()
    if ext.startswith("|") and ext.endswith("|"):
        ext = ext[1:-1]
    return {t.strip() for t in ext.split(",") if t.strip().startswith(("wU:", "wD:"))}


# ---------------------------------------------------------------------------
# data records
# ---------------------------------------------------------------------------
@dataclass
class StereoGroupRecord:
    group_type: str                  # AND | OR | ABS
    read_id: int
    atoms: List[int] = field(default_factory=list)
    bonds: List[int] = field(default_factory=list)


@dataclass
class AxialRecord:
    bond_idx: int
    chirality: str                   # cw | ccw
    ref1: int                        # neighbour of the begin atom
    ref2: int                        # neighbour of the end atom


@dataclass
class StereoInfo:
    isotopes: Dict[int, int] = field(default_factory=dict)
    abs_suffix: Dict[int, str] = field(default_factory=dict)
    local_par: Dict[int, str] = field(default_factory=dict)
    groups: List[StereoGroupRecord] = field(default_factory=list)
    axial: List[AxialRecord] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def is_empty(self) -> bool:
        return not (self.isotopes or self.abs_suffix or self.local_par
                    or self.groups or self.axial)


@dataclass
class ParsedStereoLayer:
    profile: Optional[str] = None
    cip_impl: Optional[str] = None
    tet: List[Tuple[str, str]] = field(default_factory=list)      # (atom doc id, par)
    groups: List[Tuple[str, int, List[str]]] = field(default_factory=list)
    axial: List[Tuple[str, str, str, str, str]] = field(default_factory=list)
    isotopes: List[Tuple[int, str]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def is_empty(self) -> bool:
        return not (self.tet or self.groups or self.axial or self.isotopes)


# ---------------------------------------------------------------------------
# extraction (RDKit mol -> StereoInfo)
# ---------------------------------------------------------------------------
def extract_stereo_info(
    mol: Chem.Mol,
    resolve_axial: bool = True,
    cip: Optional[Dict[int, str]] = None,
) -> StereoInfo:
    """Extract every stereo fact of a prepared mol.

    mol should come from prepare_mol (legacy perception); the modern CIP
    descriptor cache is supplied via cip or computed on a throw-away copy.
    """
    info = StereoInfo()

    # -- isotopes
    for atom in mol.GetAtoms():
        iso = atom.GetIsotope()
        if iso:
            info.isotopes[atom.GetIdx()] = int(iso)

    # -- enhanced stereo groups
    group_of_atom: Dict[int, str] = {}
    group_of_bond: Dict[int, str] = {}
    for gi, group in enumerate(mol.GetStereoGroups()):
        gname = group_type_name(group.GetGroupType())
        atoms = [a.GetIdx() for a in group.GetAtoms()]
        try:
            bonds = [b.GetIdx() for b in group.GetBonds()]
        except Exception:
            bonds = []
        info.groups.append(
            StereoGroupRecord(group_type=gname, read_id=gi + 1, atoms=atoms, bonds=bonds)
        )
        for idx in atoms:
            if idx in group_of_atom:
                info.warnings.append(
                    f"atom {idx} appears in more than one stereo group "
                    f"({group_of_atom[idx]} and {gname})"
                )
            group_of_atom[idx] = gname
        for bidx in bonds:
            group_of_bond[bidx] = gname

    # -- tetrahedral centres: local parity (source truth) + CIP suffix (cache)
    modern = cip if cip is not None else modern_cip_map(mol)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        par = tag_to_par(atom.GetChiralTag())
        if par is None:
            continue
        gname = group_of_atom.get(idx)
        cip = modern.get(idx)
        if gname in ("AND", "OR"):
            # relative configuration: parity is the only valid encoding; an
            # absolute CIP suffix here would be an over-claim (invariant I1)
            info.local_par[idx] = par
        else:
            if cip in ("R", "S", "r", "s"):
                info.abs_suffix[idx] = cip
                info.local_par[idx] = par      # V2: parity also kept
            else:
                info.warnings.append(
                    f"atom {idx} has a chiral tag but no CIP label; tag dropped"
                )

    # -- axial chirality (atropisomers)
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo not in _ATROP:
            continue
        bidx = bond.GetIdx()
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        gname = group_of_bond.get(bidx)
        if gname in ("AND", "OR"):
            info.warnings.append(
                f"axial bond {bidx} is in an {gname} stereo group; "
                "bond-level enhanced stereo is not supported yet"
            )
        stored = list(bond.GetStereoAtoms())
        ref1 = stored[0] if len(stored) >= 2 else None
        ref2 = stored[1] if len(stored) >= 2 else None
        nbrs1 = _neighbors_excluding(mol, a1, a2)
        nbrs2 = _neighbors_excluding(mol, a2, a1)
        if ref1 not in nbrs1:
            ref1 = nbrs1[0] if nbrs1 else None
        if ref2 not in nbrs2:
            ref2 = nbrs2[0] if nbrs2 else None
        if ref1 is None or ref2 is None:
            info.warnings.append(f"axial bond {bidx} has no valid reference neighbours")
            continue
        chir = "cw" if stereo == Chem.BondStereo.STEREOATROPCW else "ccw"
        info.axial.append(AxialRecord(bidx, chir, ref1, ref2))

    if resolve_axial and info.axial:
        _resolve_axial_records(mol, info)
    return info


def _mol_with_axial(mol: Chem.Mol, records: Sequence[AxialRecord]) -> Chem.Mol:
    cand = Chem.Mol(mol)
    cw, ccw = _atrop_cw_ccw()
    if cw is None or ccw is None:  # pragma: no cover - old RDKit
        return cand
    for rec in records:
        bond = cand.GetBondWithIdx(rec.bond_idx)
        bond.SetStereoAtoms(rec.ref1, rec.ref2)
        bond.SetStereo(cw if rec.chirality == "cw" else ccw)
    try:
        cand.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return cand


def _resolve_axial_records(mol: Chem.Mol, info: StereoInfo) -> None:
    """Choose (ref1, ref2, cw/ccw) that reproduce the mol's CXSMILES axis tokens.

    RDKit does not persist the reference atoms of an atropisomer, so a floating
    STEREOATROPCW cannot be written down unambiguously without choosing a
    reference frame.  We pick the frame that reproduces the molecule's canonical
    CXSMILES token, which makes the round-trip verifiable by a third-party judge.
    """
    target = _cx_axial_tokens(mol)
    if not target:
        return
    if len(info.axial) == 1:
        rec = info.axial[0]
        bond = mol.GetBondWithIdx(rec.bond_idx)
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        nbrs1 = _neighbors_excluding(mol, a1, a2)
        nbrs2 = _neighbors_excluding(mol, a2, a1)
        preferred = (rec.chirality, rec.ref1, rec.ref2)
        candidates = []
        for chir in ("cw", "ccw"):
            for r1 in nbrs1:
                for r2 in nbrs2:
                    candidates.append((chir, r1, r2))
        candidates.sort(key=lambda c: (c != preferred, c))
        for chir, r1, r2 in candidates:
            cand = _mol_with_axial(mol, [AxialRecord(rec.bond_idx, chir, r1, r2)])
            if _cx_axial_tokens(cand) == target:
                rec.chirality, rec.ref1, rec.ref2 = chir, r1, r2
                return
        info.warnings.append(
            "could not resolve axial reference frame for CXSMILES round-trip"
        )
        return

    # multiple axes: verify the extracted frames jointly, then try single flips
    if _cx_axial_tokens(_mol_with_axial(mol, info.axial)) != target:
        for i, rec in enumerate(info.axial):
            flipped = AxialRecord(
                rec.bond_idx, "ccw" if rec.chirality == "cw" else "cw",
                rec.ref1, rec.ref2,
            )
            trial = list(info.axial)
            trial[i] = flipped
            if _cx_axial_tokens(_mol_with_axial(mol, trial)) == target:
                info.axial = trial
                return
        info.warnings.append(
            "could not resolve multiple axial reference frames for CXSMILES round-trip"
        )


# ---------------------------------------------------------------------------
# rendering (StereoInfo -> MoleCode extension lines)
# ---------------------------------------------------------------------------
@dataclass
class _ObjectIdSource:
    prefix: str
    counts: Dict[str, int] = field(default_factory=dict)

    def next(self, kind: str) -> str:
        self.counts[kind] = self.counts.get(kind, 0) + 1
        return f"{self.prefix}_{kind}_{self.counts[kind]}"


def _doc_id(id_map: Dict[int, str], idx: int) -> str:
    try:
        return id_map[idx]
    except KeyError as exc:  # pragma: no cover - internal error
        raise KeyError(f"no MoleCode id for atom idx {idx}") from exc


def render_stereo_layer(
    mol: Chem.Mol,
    info: StereoInfo,
    id_map: Dict[int, str],
    prefix: str,
    parity_mode: str = "group_only",
) -> Tuple[List[str], List[str]]:
    """Render the extension layer as (header_lines, body_lines).

    Body lines are unindented.  parity_mode: all (V2, local parity for every
    centre), group_only (V1, parity only where an absolute descriptor is
    impossible), off (no S1 records).
    """
    if parity_mode not in PARITY_MODES:
        raise ValueError(f"parity_mode must be one of {PARITY_MODES}")
    header: List[str] = []
    body: List[str] = []
    if info.is_empty:
        return header, body
    ids = _ObjectIdSource(prefix)

    # S1: local parity records.  In group_only mode the absolute CIP cache is
    # enough for uppercase R/S centres (an absolute descriptor is a function of
    # the local topology), but pseudo-asymmetric r/s are *relational* (they only
    # exist when the two sides of the centre are enantiomorphic), so they need a
    # local parity record as well -- exactly like AND/OR members.
    if parity_mode != "off":
        for idx in sorted(info.local_par):
            if (
                parity_mode == "group_only"
                and info.abs_suffix.get(idx) in ("R", "S")
            ):
                continue
            body.append(
                f"{ids.next('TET')}{{tet={_doc_id(id_map, idx)}; "
                f"par={info.local_par[idx]}}}"
            )

    # S2a: enhanced stereo groups
    for rec in info.groups:
        members = ",".join(_doc_id(id_map, idx) for idx in rec.atoms)
        if not members:
            continue
        body.append(
            f"{ids.next('SG')}{{sg={rec.group_type}:{rec.read_id}; mem={members}}}"
        )

    # S2b: axial chirality
    if info.axial:
        for line in render_axial_lines(
            mol, info, id_map, prefix,
            first_index=ids.counts.get("AX", 0) + 1,
        ):
            body.append(line)
        ids.counts["AX"] = ids.counts.get("AX", 0) + len(info.axial)

    # S2c: isotopes
    for idx in sorted(info.isotopes):
        body.append(
            f"{ids.next('ISO')}{{iso={info.isotopes[idx]}; atom={_doc_id(id_map, idx)}}}"
        )

    if body:
        header = list(PROFILE_HEADER)
    return header, body


def render_axial_lines(
    mol: Chem.Mol,
    info: StereoInfo,
    id_map: Dict[int, str],
    prefix: str,
    first_index: int = 1,
) -> List[str]:
    lines: List[str] = []
    for n, rec in enumerate(info.axial, start=first_index):
        bond = mol.GetBondWithIdx(rec.bond_idx)
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        lines.append(
            f"{prefix}_AX_{n}{{ax={rec.chirality}; axis={_doc_id(id_map, a1)},"
            f"{_doc_id(id_map, a2)}; ref1={_doc_id(id_map, rec.ref1)}; "
            f"ref2={_doc_id(id_map, rec.ref2)}}}"
        )
    return lines


# ---------------------------------------------------------------------------
# parsing (MoleCode text -> ParsedStereoLayer)
# ---------------------------------------------------------------------------
def base_atom_id(doc_id: str) -> str:
    """Strip an absolute/pseudo-asymmetric CIP suffix from an atom doc id."""
    match = _SUFFIX_RE.match(doc_id.strip())
    return match.group("base") if match else doc_id.strip()


def _parse_object_body(body: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    for part in body.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        fields[key.strip()] = value.strip()
    return fields


def parse_stereo_layer(text: str) -> ParsedStereoLayer:
    """Extract the stereo/1 extension lines from a MoleCode document."""
    layer = ParsedStereoLayer()
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("%%"):
            low = line.lower()
            if "molecode-profile" in low and ":" in low:
                layer.profile = line.split(":", 1)[1].strip()
            elif "cip-impl" in low and ":" in low:
                layer.cip_impl = line.split(":", 1)[1].strip()
            continue
        match = _OBJECT_RE.match(line)
        if not match:
            continue
        oid = match.group("oid")
        fields = _parse_object_body(match.group("body"))
        try:
            if "tet" in fields and "par" in fields:
                layer.tet.append((fields["tet"], fields["par"]))
            elif "sg" in fields:
                spec = fields["sg"]
                if ":" in spec:
                    gtype, read_id = spec.split(":", 1)
                    read_id_i = int(read_id)
                else:
                    gtype, read_id_i = spec, 1
                members = [m.strip() for m in fields.get("mem", "").split(",") if m.strip()]
                layer.groups.append((gtype.strip().upper(), read_id_i, members))
            elif "ax" in fields:
                axis = [a.strip() for a in fields.get("axis", "").split(",") if a.strip()]
                if len(axis) != 2:
                    raise ValueError("ax object needs exactly two axis atom ids")
                layer.axial.append((
                    fields["ax"].strip().lower(), axis[0], axis[1],
                    fields.get("ref1", "").strip(), fields.get("ref2", "").strip(),
                ))
            elif "iso" in fields:
                mass = int(fields["iso"])
                layer.isotopes.append((mass, fields.get("atom", "").strip()))
            else:
                layer.warnings.append(f"unrecognised stereo object {oid!r}: {fields}")
        except Exception as exc:
            layer.warnings.append(f"cannot parse stereo object {oid!r}: {exc}")
    return layer


# ---------------------------------------------------------------------------
# applying the layer (ParsedStereoLayer -> RDKit mol)
# ---------------------------------------------------------------------------
def _set_chirality_from_suffix(
    mol: Chem.Mol,
    atom_idx: int,
    desired: str,
    prefer_modern: bool = True,
) -> bool:
    """Brute-force the tetrahedral tag that yields CIP label desired."""
    atom = mol.GetAtomWithIdx(atom_idx)
    for use_modern in ((True, False) if prefer_modern else (False, True)):
        for tag in (_CW, _CCW):
            if atom.HasProp("_CIPCode"):
                atom.ClearProp("_CIPCode")
            atom.SetChiralTag(tag)
            # cleanIt must stay False: a cleanIt pass removes chiral tags of
            # pseudo-asymmetric centres (legacy perception does not recognise
            # them), which would make the lowercase r/s suffix unreachable.
            try:
                Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
            except Exception:
                continue
            if use_modern and _rdCIPLabeler is not None:
                try:
                    _rdCIPLabeler.AssignCIPLabels(mol)
                except Exception:
                    pass
            got: Optional[str] = None
            if atom.HasProp("_CIPCode"):
                got = atom.GetProp("_CIPCode")
            if got is None:
                continue
            if got == desired or (not use_modern and got.upper() == desired.upper()):
                return True
    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return False


def apply_stereo_layer(
    mol: Chem.Mol,
    layer: ParsedStereoLayer,
    doc_id_to_idx: Dict[str, int],
) -> Tuple[Chem.Mol, List[str]]:
    """Apply a parsed extension layer to a base-parsed mol.

    Returns (new_mol, issues) where issues is a list of human-readable strings
    (empty when the document is consistent).
    """
    issues: List[str] = list(layer.warnings)
    rw = Chem.RWMol(mol)

    def resolve(doc_id: str, what: str = "atom") -> Optional[int]:
        idx = doc_id_to_idx.get(doc_id)
        if idx is None:
            issues.append(f"unknown {what} {doc_id!r} in stereo layer")
        return idx

    # I4: isotopes first (CIP labels depend on them)
    for mass, doc_id in layer.isotopes:
        idx = resolve(doc_id, "isotope atom")
        if idx is None:
            continue
        rw.GetAtomWithIdx(idx).SetIsotope(mass)

    # S1: local parity records (source of truth)
    for doc_id, par in layer.tet:
        idx = resolve(doc_id)
        if idx is None:
            continue
        try:
            tag = par_to_tag(par)
        except ValueError as exc:
            issues.append(str(exc))
            continue
        rw.GetAtomWithIdx(idx).SetChiralTag(tag)

    # S2a: enhanced stereo groups.  The group objects must be created on (and
    # attached to) the same RWMol instance, otherwise RDKit silently drops them.
    groups = []
    group_members: Dict[int, str] = {}
    for gtype, read_id, member_ids in layer.groups:
        if gtype not in GROUP_NAME_TO_TYPE:
            issues.append(f"unknown stereo group type {gtype!r}")
            continue
        idxs: List[int] = []
        for doc_id in member_ids:
            idx = resolve(doc_id)
            if idx is None:
                continue
            if idx in group_members:
                issues.append(
                    f"atom {doc_id} belongs to two stereo groups "
                    f"({group_members[idx]} and {gtype})"
                )
            group_members[idx] = gtype
            idxs.append(idx)
        if not idxs:
            continue
        if not hasattr(Chem, "CreateStereoGroup"):  # pragma: no cover - old RDKit
            issues.append(
                "this RDKit build has no Chem.CreateStereoGroup; "
                "enhanced-stereo groups are kept as local parity only"
            )
            groups = []
            break
        try:
            groups.append(Chem.CreateStereoGroup(GROUP_NAME_TO_TYPE[gtype], rw, idxs, []))
        except Exception as exc:
            issues.append(f"cannot build stereo group {gtype}:{read_id}: {exc}")
    if groups:
        try:
            rw.SetStereoGroups(groups)
        except Exception as exc:
            issues.append(f"cannot attach stereo groups: {exc}")

    # S2b: axial chirality
    for chir, id1, id2, ref1, ref2 in layer.axial:
        i1, i2 = resolve(id1), resolve(id2)
        r1, r2 = resolve(ref1, "axial ref1"), resolve(ref2, "axial ref2")
        if None in (i1, i2, r1, r2):
            continue
        bond = rw.GetBondBetweenAtoms(i1, i2)
        if bond is None:
            issues.append(f"axial axis {id1}-{id2} is not a bond")
            continue
        nbr1 = {
            n.GetIdx()
            for n in rw.GetAtomWithIdx(i1).GetNeighbors()
            if n.GetIdx() != i2
        }
        nbr2 = {
            n.GetIdx()
            for n in rw.GetAtomWithIdx(i2).GetNeighbors()
            if n.GetIdx() != i1
        }
        if r1 not in nbr1:
            issues.append(f"axial ref1 {ref1} is not a neighbour of {id1}")
            continue
        if r2 not in nbr2:
            issues.append(f"axial ref2 {ref2} is not a neighbour of {id2}")
            continue
        cw, ccw = _atrop_cw_ccw()
        if cw is None or ccw is None:  # pragma: no cover - old RDKit
            issues.append(
                f"axial stereo on {id1}-{id2} requires an RDKit build with "
                "STEREOATROPCW/STEREOATROPCCW; record ignored"
            )
            continue
        try:
            bond.SetStereoAtoms(r1, r2)
            bond.SetStereo(cw if chir == "cw" else ccw)
        except Exception as exc:
            issues.append(f"cannot set axial stereo on {id1}-{id2}: {exc}")

    mol1 = rw.GetMol()

    # I4 continuation: CIP labels re-derived on the final graph
    if _rdCIPLabeler is not None:
        try:
            _rdCIPLabeler.AssignCIPLabels(mol1)
        except Exception:
            pass
    try:
        Chem.AssignStereochemistry(mol1, cleanIt=False, force=True)
    except Exception:
        pass

    return mol1, issues


# ---------------------------------------------------------------------------
# fingerprints and validation
# ---------------------------------------------------------------------------
def _group_signature(mol: Chem.Mol) -> List[Tuple[str, Tuple[int, ...]]]:
    sig = []
    for group in mol.GetStereoGroups():
        sig.append((
            group_type_name(group.GetGroupType()),
            tuple(sorted(a.GetIdx() for a in group.GetAtoms())),
        ))
    return sorted(sig)


def _axial_signature(mol: Chem.Mol) -> List[Tuple[int, int, str]]:
    out = []
    for bond in mol.GetBonds():
        if bond.GetStereo() in _ATROP:
            out.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), str(bond.GetStereo())))
    return sorted(out)


def _isotope_signature(mol: Chem.Mol) -> List[Tuple[int, int]]:
    return sorted((a.GetIdx(), a.GetIsotope()) for a in mol.GetAtoms() if a.GetIsotope())


def stereo_fingerprint(mol: Chem.Mol) -> Dict[str, object]:
    """Index-based stereo fingerprint (for tests of the same atom ordering)."""
    m = prepare_mol(mol)
    return {
        "smiles": canonical_smiles(m),
        "inchi_key": inchi_key(m),
        "cxsmiles": cx_smiles(m),
        "cip": dict(sorted(modern_cip_map(mol).items())),
        "groups": _group_signature(m),
        "axial": _axial_signature(m),
        "isotopes": _isotope_signature(m),
    }


def validate_document(
    text: str,
    suffix_map: Optional[Dict[str, str]] = None,
) -> Dict[str, object]:
    """Static checks on a MoleCode stereo document (no RDKit required).

    suffix_map maps atom doc id -> CIP suffix as parsed from the atom ids; it is
    passed in by the parser (which strips the suffix when building the mol).
    Returns a dict with ok and an issues list.
    """
    issues: List[str] = []
    layer = parse_stereo_layer(text)
    suffix_map = suffix_map or {}

    # I2: duplicate / invalid local parity records
    seen_tet: Dict[str, str] = {}
    for doc_id, par in layer.tet:
        if par.strip().lower() not in ("cw", "ccw"):
            issues.append(f"I2 invalid local parity {par!r} for atom {doc_id}")
        if doc_id in seen_tet:
            issues.append(
                f"I2 duplicate TET record for atom {doc_id} "
                f"({seen_tet[doc_id]} then {par}); last one would win"
            )
        seen_tet[doc_id] = par

    # I3: each stereo element belongs to at most one group
    membership: Dict[str, List[str]] = {}
    for gtype, _read_id, members in layer.groups:
        if gtype not in ("AND", "OR", "ABS"):
            issues.append(f"I3 unknown stereo group type {gtype!r}")
        for member in members:
            membership.setdefault(member, []).append(gtype)
    for atom_id, gtypes in sorted(membership.items()):
        if len(gtypes) > 1:
            issues.append(f"I3 group conflict: atom {atom_id} is in groups {gtypes}")

    # I1: no absolute suffix on AND/OR members (match by base id so that a
    # tampered "C_1_S" is caught even when the group member is written "C_1")
    group_of_atom: Dict[str, str] = {}
    for gtype, _read_id, members in layer.groups:
        for member in members:
            group_of_atom[member] = gtype
            group_of_atom[base_atom_id(member)] = gtype
    for atom_id, cip in suffix_map.items():
        gtype = group_of_atom.get(atom_id) or group_of_atom.get(base_atom_id(atom_id))
        if gtype in ("AND", "OR"):
            issues.append(
                f"I1 over-claim: atom {atom_id} has absolute suffix _{cip} "
                f"but belongs to an {gtype} group"
            )

    # I3: ABS and AND/OR membership must be disjoint
    abs_members = {m for g, _r, ms in layer.groups if g == "ABS" for m in ms}
    rel_members = {m for g, _r, ms in layer.groups if g in ("AND", "OR") for m in ms}
    for atom_id in sorted(abs_members & rel_members):
        issues.append(f"I3 group conflict: atom {atom_id} is both ABS and relative")

    # I5: axial refs must be present
    for chir, id1, id2, ref1, ref2 in layer.axial:
        if not ref1 or not ref2:
            issues.append(f"I5 axial {id1}-{id2} is missing ref1/ref2")
        if chir not in ("cw", "ccw"):
            issues.append(f"I5 axial {id1}-{id2} has invalid chirality {chir!r}")

    # I6: judge adequacy
    if any(g in ("AND", "OR") for g, _r, _m in layer.groups) or layer.axial:
        issues.append(
            "I6 judge: full InChIKey is not a valid judge here "
            "(InChI drops AND/OR groups and is blind to axial chirality)"
        )

    hard = [i for i in issues if not i.startswith("I6")]
    return {"ok": not hard, "issues": issues,
            "profile": layer.profile, "cip_impl": layer.cip_impl}


# Public alias (used by the reader).  Kept at the end so the private helper can
# stay close to apply_stereo_layer.
set_chirality_from_suffix = _set_chirality_from_suffix


def stereo_coverage(mol: Chem.Mol) -> Dict[str, int]:
    """Count stereo elements of a molecule (used by the benchmark harness)."""
    m = prepare_mol(mol)
    cip = cip_map(m)
    return {
        "centres": len(cip),
        "pseudo_asymmetric": sum(1 for v in cip.values() if v in ("r", "s")),
        "groups": len(m.GetStereoGroups()),
        "axial": sum(1 for b in m.GetBonds() if b.GetStereo() in _ATROP),
        "isotopes": sum(1 for a in m.GetAtoms() if a.GetIsotope()),
    }

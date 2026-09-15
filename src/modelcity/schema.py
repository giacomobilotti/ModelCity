"""Column-role resolution and validation.

Datasets name their columns differently: Yautepec calls the identifier ``Sitio``
and the size measure ``Area``; Viabundus calls them ``Nodes_ID`` and
``Inhabitants``.  Rather than force one naming convention, the caller maps
*roles* to whatever the columns happen to be called.

Resolution is deliberately predictable: a role is satisfied either by an
explicit mapping or by a column already carrying the canonical role name.
Anything else is reported as an error with a suggestion, never guessed at.
"""

from __future__ import annotations

import difflib
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = [
    "SchemaError",
    "ColumnRoles",
    "REQUIRED_ROLES",
    "OPTIONAL_ROLES",
    "ALL_ROLES",
    "resolve_roles",
    "suggest_column",
]


REQUIRED_ROLES: Tuple[str, ...] = ("id", "size", "start")
OPTIONAL_ROLES: Tuple[str, ...] = ("end", "name", "group")
ALL_ROLES: Tuple[str, ...] = REQUIRED_ROLES + OPTIONAL_ROLES

ROLE_HELP: Dict[str, str] = {
    "id": "identifier of the settlement, repeated across periods",
    "size": "numeric size measure (area or population)",
    "start": "start of the period, or the observation year for snapshot data",
    "end": "end of the period; omit for snapshot data",
    "name": "human-readable label for the settlement",
    "group": "grouping such as ceramic phase or region",
}

# Used only to suggest columns in error messages.  Never used to resolve a role
# silently: a wrong guess about which column holds the size measure would
# quietly corrupt every downstream figure.
ROLE_ALIASES: Dict[str, Tuple[str, ...]] = {
    "id": (
        "id", "ids", "site", "site_id", "siteid", "sitio", "sitios",
        "node", "node_id", "nodes_id", "city_id", "settlement", "settlement_id",
    ),
    "size": (
        "size", "area", "areas", "inhabitants", "population", "pop",
        "hectares", "ha", "extent", "surface", "people",
    ),
    "start": (
        "start", "start_bp", "start_year", "startyear", "from", "begin",
        "beginning", "year", "date", "period_start", "t_start",
    ),
    "end": (
        "end", "end_bp", "end_year", "endyear", "to", "stop", "finish",
        "period_end", "t_end",
    ),
    "name": ("name", "names", "label", "toponym", "site_name", "city", "place"),
    "group": (
        "group", "groups", "region", "regions", "fase", "phase", "period",
        "culture", "category",
    ),
}


class SchemaError(ValueError):
    """Raised when the input columns cannot be mapped onto the canonical roles."""


@dataclass(frozen=True)
class ColumnRoles:
    """The resolved mapping from canonical role to actual column name."""

    id: str
    size: str
    start: str
    end: Optional[str] = None
    name: Optional[str] = None
    group: Optional[str] = None

    def mapping(self) -> Dict[str, str]:
        """Role to column name, omitting roles that were not supplied."""
        pairs = (
            ("id", self.id),
            ("size", self.size),
            ("start", self.start),
            ("end", self.end),
            ("name", self.name),
            ("group", self.group),
        )
        return {role: column for role, column in pairs if column is not None}

    @property
    def has_end(self) -> bool:
        return self.end is not None


def _tokenise(name: str) -> List[str]:
    """Split a column name into words, so ``start_bp`` yields ``start`` and ``bp``."""
    return [token for token in re.split(r"[^a-z0-9]+", name.lower()) if token]


def suggest_column(role: str, columns: Sequence[str]) -> Optional[str]:
    """Best guess at which column fills ``role``, or ``None`` if nothing fits.

    Tries the role's known aliases first, since ``size`` and ``Area`` share no
    characters and so are invisible to pure string similarity.
    """
    if not columns:
        return None

    lowered = {str(column).lower(): str(column) for column in columns}
    aliases = ROLE_ALIASES.get(role, (role,))

    for alias in aliases:
        if alias in lowered:
            return lowered[alias]

    # Fall back to fuzzy matching, against the aliases as well as the role name
    # so that e.g. "Inhabitants_1500" still points at the size role.
    for candidate_set in (aliases, (role,)):
        for alias in candidate_set:
            close = difflib.get_close_matches(alias, list(lowered), n=1, cutoff=0.7)
            if close:
                return lowered[close[0]]

    # Last resort: a compound name such as "Inhabitants_1500" that contains an
    # alias as a whole word.  Matching on bare substrings instead would let a
    # short alias like "ha" fire on any column containing those letters.
    for alias in aliases:
        for lower_name, original in lowered.items():
            if alias in _tokenise(lower_name):
                return original

    return None


def _resolve_one(
    role: str,
    supplied: Optional[str],
    columns: Sequence[str],
) -> Tuple[Optional[str], Optional[str]]:
    """Resolve a single role.

    Returns ``(resolved_column, problem)`` where exactly one is not ``None``.
    """
    column_list = [str(column) for column in columns]

    if supplied is not None:
        if supplied in column_list:
            return supplied, None
        return None, "missing"

    if role in column_list:
        return role, None

    # Accept a case-only difference, which is common when data is exported from
    # different tools, but nothing looser than that.
    for column in column_list:
        if column.lower() == role:
            return column, None

    return None, "unmapped"


def _format_problem(
    role: str,
    supplied: Optional[str],
    problem: str,
    columns: Sequence[str],
    suggestion: Optional[str],
) -> str:
    column_list = [str(column) for column in columns]
    if problem == "missing":
        headline = (
            "column {!r} was supplied for role {!r} but is not in the data.".format(
                supplied, role
            )
        )
    else:
        headline = "required column for role {!r} was not found.".format(role)

    lines = [
        headline,
        "  role means: {}".format(ROLE_HELP.get(role, "")),
        "  looked for: {!r}".format(supplied if supplied is not None else role),
        "  available:  {}".format(", ".join(column_list) if column_list else "<none>"),
    ]
    if suggestion is not None:
        lines.append("  did you mean {!r}?".format(suggestion))
    return "\n".join(lines)


def _format_fix(resolved: Dict[str, str], guesses: Dict[str, Optional[str]], size_type: Optional[str]) -> str:
    """Build a copy-pasteable ``as_settlements`` call.

    Roles that resolved keep their real column; roles that failed use the
    suggestion, or a placeholder when there is nothing sensible to suggest.
    """
    args: List[str] = []
    for role in ALL_ROLES:
        if role in resolved:
            args.append('{}="{}"'.format(role, resolved[role]))
        elif role in guesses:
            guess = guesses[role]
            args.append('{}="{}"'.format(role, guess if guess else "<column>"))
    args.append('size_type="{}"'.format(size_type if size_type else "<area|population>"))
    return "  fix: as_settlements(df, {})".format(", ".join(args))


def resolve_roles(
    columns: Sequence[str],
    *,
    id: Optional[str] = None,
    size: Optional[str] = None,
    start: Optional[str] = None,
    end: Optional[str] = None,
    name: Optional[str] = None,
    group: Optional[str] = None,
    size_type: Optional[str] = None,
) -> ColumnRoles:
    """Map roles onto ``columns``, raising :class:`SchemaError` if any required
    role cannot be satisfied.

    All failures are collected into a single exception rather than surfacing one
    at a time, so a mislabelled dataset can be fixed in one pass.
    """
    column_list = [str(column) for column in columns]
    supplied = {
        "id": id,
        "size": size,
        "start": start,
        "end": end,
        "name": name,
        "group": group,
    }

    resolved: Dict[str, str] = {}
    problems: List[str] = []
    guesses: Dict[str, Optional[str]] = {}

    for role in ALL_ROLES:
        value = supplied[role]
        column, problem = _resolve_one(role, value, column_list)

        if column is not None:
            resolved[role] = column
            continue

        # An optional role that was never requested is simply absent.
        if role in OPTIONAL_ROLES and value is None:
            continue

        suggestion = suggest_column(role, column_list)
        guesses[role] = suggestion
        problems.append(_format_problem(role, value, problem, column_list, suggestion))

    if problems:
        body = "\n".join(problems)
        raise SchemaError(
            "{}\n{}".format(body, _format_fix(resolved, guesses, size_type))
        )

    return ColumnRoles(
        id=resolved["id"],
        size=resolved["size"],
        start=resolved["start"],
        end=resolved.get("end"),
        name=resolved.get("name"),
        group=resolved.get("group"),
    )

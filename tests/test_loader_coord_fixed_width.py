"""Regression: AMBER restart coordinates must be parsed FIXED-WIDTH (F12.7),
not whitespace-split.

Reviewer finding (fast-amber-loader rewrite). ``amber_loader.
read_amber_coordinates`` reads the coordinate block with ``coord_text.split()``
(whitespace), but the AMBER restart/inpcrd format is a fixed-width Fortran
``6F12.7`` record with NO guaranteed separator between fields. Any coordinate
that fills all 12 columns -- e.g. a negative value with >= 3 integer digits
(``"%12.7f" % -100.0 == "-100.0000000"``, exactly 12 chars) or a positive
value with >= 4 integer digits (``"1000.0000000"``) -- abuts its neighbour with
zero intervening whitespace. ``split()`` then fuses adjacent fields into one
unparseable token and the read raises, even though OpenMM's own
``AmberAsciiRestart`` (the module ``amber_loader`` documents itself as
mirroring) and ParmEd both parse the file correctly by column slicing.

Such coordinate magnitudes are entirely normal for the large, unwrapped
systems this project targets (up to 1M atoms). The failure is fail-loud (a
``ValueError``, not silent corruption), but it rejects well-formed input that
the reference readers accept -- a correctness regression versus the loader's
stated reference.

Expected fix: parse the coordinate block by 12-column slices (the module
already does exactly this for the box line via ``_read_fixed_width_floats``);
this test then passes.
"""

from __future__ import annotations

import pytest

pytest.importorskip("robosample")

from robosample import amber_loader  # noqa: E402


def _write_restart(path, natom: int, value: float) -> None:
    field = "%12.7f" % value
    assert len(field) == 12, f"test setup: field {field!r} is not 12 wide"
    n_coord_lines = (natom + 1) // 2
    with open(path, "w") as fh:
        fh.write("packed-coordinate restart\n")
        fh.write("%5d\n" % natom)
        remaining = natom
        for _ in range(n_coord_lines):
            per_line = min(2, remaining)
            fh.write(field * (3 * per_line) + "\n")
            remaining -= per_line


def test_packed_fixed_width_coordinates(tmp_path):
    """A restart whose F12.7 fields abut with no whitespace must still parse."""
    path = tmp_path / "packed.rst7"
    # 4 atoms => 2 coordinate lines; -123.456789 formats to exactly 12 chars,
    # so all six fields per line run together with zero spaces.
    _write_restart(path, natom=4, value=-123.456789)

    coords, box = amber_loader.read_amber_coordinates(str(path))

    assert box is None
    assert coords.shape == (4, 3)
    # -123.456789 A -> nm (ANG_TO_NM = 0.1)
    expected = -123.456789 * amber_loader.ANG_TO_NM
    assert coords == pytest.approx(expected)

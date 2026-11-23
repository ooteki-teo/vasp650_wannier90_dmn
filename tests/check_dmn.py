"""Check a VASP DMN against the input crystal, EIG, AMN, and optionally QE."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
import spglib


def read_dmn(path, nw):
    with path.open("rb") as stream:
        stream.readline()
        nb, ns, ni, nk = map(int, stream.readline().split())
        head = bytearray()
        while b"(" not in (line := stream.readline()):
            if not line:
                raise ValueError(f"No complex matrices in {path}")
            head.extend(line)
        maps = np.fromstring(head.decode(), sep=" ", dtype=int) - 1
        assert len(maps) == nk + ni + ni * ns
        raw = (line + stream.read()).translate(bytes.maketrans(b"(,)", b"   "))
        data = np.fromstring(raw.decode(), sep=" ").view(np.complex128)
        assert data.size == ni * ns * (nw * nw + nb * nb)
    offset = ni * ns * nw * nw
    dw = data[:offset].reshape(ni, ns, nw, nw).swapaxes(-1, -2)
    db = data[offset:].reshape(ni, ns, nb, nb).swapaxes(-1, -2)
    return maps[:nk], maps[nk:nk + ni], maps[nk + ni:].reshape(ni, ns), dw, db


def read_win(path):
    text = path.read_text()
    kblock = re.search(r"(?ims)^\s*begin kpoints\s*\n(.*?)^\s*end kpoints", text)
    nw = int(re.search(r"(?im)^\s*num_wann\s*=\s*(\d+)", text)[1])
    return np.loadtxt(kblock[1].splitlines()), nw


def read_operations(path):
    fields = [s.split()[2:] for s in path.read_text().splitlines()
              if s.startswith(" DMN operation")]
    assert fields
    nops = sum(int(row[0]) == 1 for row in fields)
    fields = fields[:len(fields) // nops]
    rotation = np.array([row[1:10] for row in fields], dtype=int).reshape(-1, 3, 3)
    translation = np.array([row[10:13] for row in fields], dtype=float)
    return rotation, translation


def read_structure(path):
    lines = path.read_text().splitlines()
    lattice = np.loadtxt(lines[2:5]) * float(lines[1])
    counts = list(map(int, lines[6].split()))
    labels = np.repeat(np.arange(1, len(counts) + 1), counts)
    positions = np.loadtxt(lines[8:8 + sum(counts)], usecols=(0, 1, 2))
    if lines[7].lower().startswith("c"):
        positions = positions * float(lines[1]) @ np.linalg.inv(lattice)
    return lattice, positions, labels


def wrap(vector):
    return vector - np.rint(vector)


def signature(kpoints):
    return sorted(tuple(np.round(k % 1, 8) % 1) for k in kpoints)


def check(seed, magnetic, qe_seed, window):
    path = seed.parent
    base_win = path / "wannier90.win"
    k, nw = read_win(base_win)
    data = read_dmn(seed.with_suffix(seed.suffix + ".dmn"), nw)
    k_to_ir, ir_to_k, images, dw, db = data
    ni, ns = images.shape
    nb = db.shape[-1]
    rotations, translations = read_operations(path / "OUTCAR")
    assert len(rotations) == ns
    cell = read_structure(path / "POSCAR")
    if magnetic:
        operations = spglib.get_magnetic_symmetry((*cell, [2., -2., 0., 0.]))
        mask = ~operations["time_reversals"]
    else:
        operations = spglib.get_symmetry(cell)
        mask = np.ones(len(operations["rotations"]), dtype=bool)
    reference_r = operations["rotations"][mask]
    reference_t = operations["translations"][mask]
    assert len(reference_r) == ns
    for rotation, translation in zip(rotations, translations, strict=True):
        match_r = (reference_r == rotation).all(axis=(1, 2))
        match_t = (abs(wrap(reference_t - translation)) < 1e-6).all(axis=1)
        assert np.count_nonzero(match_r & match_t) == 1

    rk = np.rint(np.linalg.inv(rotations)).astype(int)
    mapped = np.einsum("ia,sab->isb", k[ir_to_k], rk)
    k_error = np.max(abs(wrap(mapped - k[images])))
    assert np.array_equal(k_to_ir[images], np.broadcast_to(np.arange(ni)[:, None], images.shape))
    assert np.array_equal(images[:, 0], ir_to_k)
    assert len(np.unique(images)) == len(k)
    err_dw = np.max(abs(dw.conj().swapaxes(-1, -2) @ dw - np.eye(nw)))
    err_id = np.max(abs(db[:, 0] - np.eye(nb)))

    e = np.loadtxt(str(seed) + ".eig")[:, -1].reshape(len(k), nb)
    a = np.loadtxt(str(seed) + ".amn", skiprows=2)[:, -2:].copy().view(np.complex128)
    a = a.reshape(len(k), nw, nb).swapaxes(-1, -2)
    err_e = np.max(abs((e[images][..., :, None] - e[ir_to_k, None, None, :]) * db))
    delta_a = db @ a[ir_to_k, None] - a[images] @ dw
    scale = np.linalg.norm(a[images], axis=(-1, -2))
    err_a = np.max(np.linalg.norm(delta_a, axis=(-1, -2)) / scale)
    err_a_out = np.max(abs(delta_a))
    err_norm = np.max(abs(1 - np.sum(abs(db)**2, axis=-2)))

    result = dict(seed=str(seed), nbands=nb, nwann=nw, operations=ns, irreducible=ni,
                  nkpoints=len(k), k_map_max=float(k_error), dw_unitarity_max=float(err_dw),
                  db_identity_max=float(err_id), db_band_norm_max=float(err_norm),
                  e_sewing_max_eV=float(err_e), amn_sewing_relative_max=float(err_a),
                  amn_sewing_element_max=float(err_a_out))
    if window is not None:
        selected = e <= window
        rows = selected[images]
        cols = np.broadcast_to(selected[ir_to_k, None, :], (*images.shape, nb))
        assert np.array_equal(rows.sum(axis=-1), cols.sum(axis=-1))
        db_window = db * rows[..., :, None] * cols[..., None, :]
        cropped_a = a * selected[..., None]
        delta = db_window @ cropped_a[ir_to_k, None] - cropped_a[images] @ dw
        result["window"] = dict(
            max_eV=window, min_states=int(selected.sum(axis=1).min()),
            band_norm_max=float(np.max(abs(1 - np.sum(abs(db_window)**2, axis=-2))[cols])),
            e_sewing_max_eV=float(np.max(abs((e[images][..., :, None] - e[ir_to_k, None, None, :])
                                            * db_window))),
            amn_sewing_relative_max=float(np.max(np.linalg.norm(delta, axis=(-1, -2))
                                                / np.linalg.norm(cropped_a[images], axis=(-1, -2)))))
    if qe_seed is not None:
        kq, nwq = read_win(Path(str(qe_seed) + ".win"))
        _, qrep, qi, qw, qb = read_dmn(Path(str(qe_seed) + ".dmn"), nwq)
        assert qi.shape == images.shape
        orbits = {tuple(signature(k[np.unique(img)])): i for i, img in enumerate(images)}
        for img in qi:
            assert tuple(signature(kq[np.unique(img)])) in orbits
        ig, qg = np.argmin(np.linalg.norm(k[ir_to_k], axis=1)), np.argmin(np.linalg.norm(kq[qrep], axis=1))
        err_t = {}
        for key, x, y in (("wannier", dw[ig], qw[qg]), ("bloch", db[ig], qb[qg])):
            err_t[key] = float(np.max(abs(np.sort_complex(np.trace(x, axis1=1, axis2=2).round(7))
                                        - np.sort_complex(np.trace(y, axis1=1, axis2=2).round(7)))))
        result["qe_gamma_character_max"] = err_t
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=Path)
    parser.add_argument("--magnetic", action="store_true")
    parser.add_argument("--qe-seed", type=Path)
    parser.add_argument("--window", type=float)
    parser.add_argument("--check", action="store_true", help="Fail on numerical errors")
    args = parser.parse_args()
    report = check(args.seed, args.magnetic, args.qe_seed, args.window)
    print(json.dumps(report, indent=2))
    if args.check:
        assert report["k_map_max"] < 1e-7
        assert report["dw_unitarity_max"] < 1e-7
        assert report["db_identity_max"] < 1e-7
        metrics = report.get("window", report)
        assert metrics["e_sewing_max_eV"] < 1e-5
        assert metrics["amn_sewing_relative_max"] < (5e-3 if args.magnetic else 1e-5)
        assert metrics.get("band_norm_max", report["db_band_norm_max"]) < 1e-6
        assert max(report.get("qe_gamma_character_max", {"none": 0}).values()) < 1e-4

"""Locate TaAs Weyl nodes independently near the symmetry images of W1 and W2."""
import argparse
import json
from pathlib import Path
import re

import numpy as np
from scipy.optimize import minimize

from check_dmn import read_operations, read_structure, wrap


class Hamiltonian:
    def __init__(self, seed):
        with Path(f"{seed}_hr.dat").open() as stream:
            stream.readline()
            nw = int(stream.readline())
            nr = int(stream.readline())
            degeneracy = np.fromstring("".join(stream.readline() for _ in range((nr + 14) // 15)),
                                       sep=" ", dtype=int)
            raw = np.loadtxt(stream)
        assert len(degeneracy) == nr and len(raw) == nr * nw**2
        hopping = {
            tuple(map(int, row[:5])): complex(*row[5:]) / deg
            for row, deg in zip(raw, np.repeat(degeneracy, nw**2), strict=True)
        }
        terms = {}
        with Path(f"{seed}_wsvec.dat").open() as stream:
            stream.readline()
            while line := stream.readline():
                key = tuple(map(int, line.split()))
                r = np.array(key[:3])
                i, j = key[3] - 1, key[4] - 1
                n = int(stream.readline())
                amplitude = hopping.pop(key) / n
                for _ in range(n):
                    shift = np.fromstring(stream.readline(), dtype=int, sep=" ")
                    vector = tuple(r + shift)
                    if vector not in terms:
                        terms[vector] = np.zeros((nw, nw), dtype=complex)
                    terms[vector][i, j] += amplitude
        assert not hopping
        self.r = np.asarray(list(terms))
        self.h = np.asarray(list(terms.values())).reshape(-1, nw**2)
        self.nw = nw

    def evaluate(self, k, derivatives=False):
        phase = np.exp(2j * np.pi * (self.r @ k))
        matrix = (phase @ self.h).reshape(self.nw, self.nw)
        if derivatives:
            dh = ((2j * np.pi * self.r.T * phase) @ self.h).reshape(3, self.nw, self.nw)
            return matrix, dh
        return matrix


def unique_images(seed, rk, lattice, tolerance=1e-5):
    images = []
    for sign in (1, -1):
        for i, rot in enumerate(rk):
            k = sign * (seed @ rot)
            distance = [np.linalg.norm(wrap(k - p) @ lattice) for p, _, _ in images]
            if not distance or min(distance) > tolerance:
                images.append((k, i, sign))
    return images


def solve_node(ham, start, occupied, lattice):
    def gap_square(k):
        h, dh = ham.evaluate(k, derivatives=True)
        e, v = np.linalg.eigh(h)
        c = v[:, occupied - 1:occupied + 1]
        slopes = np.einsum("ai,pab,bi->pi", c.conj(), dh, c).real
        gap = e[occupied] - e[occupied - 1]
        return gap**2, 2 * gap * (slopes[:, 1] - slopes[:, 0])

    result = minimize(gap_square, start, method="BFGS", jac=True,
                      options={"gtol": 2e-14, "maxiter": 250})
    k = wrap(result.x)
    matrix, dh = ham.evaluate(k, derivatives=True)
    e, v = np.linalg.eigh(matrix)
    c = v[:, occupied - 1:occupied + 1]
    dh_cart = np.einsum("pq,qab->pab", lattice.T / (2 * np.pi), dh)
    d = np.einsum("ai,pab,bj->pij", c.conj(), dh_cart, c)
    pauli = np.array([d[:, 0, 1].real, -d[:, 0, 1].imag,
                      0.5 * (d[:, 0, 0] - d[:, 1, 1]).real])
    det = np.linalg.det(pauli)
    isolation = min(e[occupied - 1] - e[occupied - 2], e[occupied + 1] - e[occupied])
    return k, dict(gap_eV=float(e[occupied] - e[occupied - 1]),
                   energy_eV=float((e[occupied] + e[occupied - 1]) / 2),
                   isolation_eV=float(isolation), velocity_det=float(det),
                   chirality=int(np.sign(det)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--fermi", type=float, help="Defaults to the final E-fermi in OUTCAR")
    args = parser.parse_args()
    ham = Hamiltonian(args.seed)
    lattice, _, _ = read_structure(args.seed.parent / "POSCAR")
    recip = 2 * np.pi * np.linalg.inv(lattice).T
    outcar = args.seed.parent / "OUTCAR"
    rotation, _ = read_operations(outcar)
    fermi = args.fermi
    if fermi is None:
        fermi = float(re.findall(r"E-fermi\s*:\s*([-+.\d]+)", outcar.read_text())[-1])
    rk = np.rint(np.linalg.inv(rotation)).astype(int)
    conventional_lengths = np.array([lattice[0, 0], lattice[1, 1], 2 * lattice[2, 2]])

    # Fractional conventional coordinates from Garcia et al., Phys. Rev. Research 2, 013073.
    seeds = {"W1": [0.008, 0.505, 0], "W2": [0.020, 0.280, 0.578]}
    output = {}
    for label, point in seeds.items():
        start = (np.array(point) / conventional_lengths) @ lattice.T
        representative, _ = solve_node(ham, start, 16, lattice)
        # Use crystallography to predict the images and independently minimize around each one.
        expected = unique_images(start + wrap(representative - start), rk, recip)
        nodes = []
        ref_chirality = None
        for k_expected, operation, time_reversal_sign in expected:
            offset = np.array([1.2e-4, -1.5e-4, 1.e-4])
            k, record = solve_node(ham, k_expected + offset, 16, lattice)
            if ref_chirality is None:
                ref_chirality = record["chirality"]
            k_unfolded = k_expected + wrap(k - k_expected)
            record.update(frac_primitive=k.tolist(),
                          frac_conventional=((k_unfolded @ recip) * conventional_lengths / (2 * np.pi)).tolist(),
                          sym_error_A_inv=float(np.linalg.norm(wrap(k - k_expected) @ recip)),
                          sym_image_gap_eV=float(np.diff(np.linalg.eigvalsh(ham.evaluate(k_expected))[15:17])[0]),
                          expected_chirality=int(ref_chirality * round(np.linalg.det(rotation[operation]))),
                          point_operation=operation + 1, time_reversal=bool(time_reversal_sign < 0))
            nodes.append(record)
        e = np.array([n["energy_eV"] for n in nodes])
        positions = np.array([n["frac_primitive"] for n in nodes])
        distances = np.linalg.norm(wrap(positions[:, None] - positions[None, :]) @ recip, axis=-1)
        nearest = distances[~np.eye(len(nodes), dtype=bool)]
        output[label] = dict(count=len(nodes), closest_pair_A_inv=float(nearest.min()) if nearest.size else None, max_gap_eV=max(n["gap_eV"] for n in nodes),
                             max_sym_error_A_inv=max(n["sym_error_A_inv"] for n in nodes),
                             max_sym_image_gap_eV=max(n["sym_image_gap_eV"] for n in nodes),
                             energy_range_eV=float(np.ptp(e)),
                             mean_energy_relative_to_fermi_eV=float(np.mean(e) - fermi),
                             nodes=nodes)
    # Sample generic points as a separate check of the interpolated spectra.
    rng = np.random.default_rng(105)
    covariance = {"unitary_max_eV": 0.0, "time_reversal_max_eV": 0.0}
    for k in rng.uniform(-0.5, 0.5, size=(12, 3)):
        e = np.linalg.eigvalsh(ham.evaluate(k))
        covariance["time_reversal_max_eV"] = max(covariance["time_reversal_max_eV"],
                                               float(np.max(abs(e - np.linalg.eigvalsh(ham.evaluate(-k))))))
        for op in rk:
            delta = np.linalg.eigvalsh(ham.evaluate(k @ op)) - e
            covariance["unitary_max_eV"] = max(covariance["unitary_max_eV"], float(np.max(abs(delta))))
    output["spectra"] = covariance
    print(json.dumps(output, indent=2))
    if args.check:
        kz = np.array([n["frac_conventional"][2] for n in output["W1"]["nodes"]])
        assert max(abs(wrap(kz))) < 1e-5
        for label, expected_count in (("W1", 8), ("W2", 16)):
            orbit = output[label]
            assert orbit["count"] == expected_count
            assert orbit["closest_pair_A_inv"] > 1e-4
            assert orbit["max_gap_eV"] < 1e-6
            assert orbit["max_sym_image_gap_eV"] < 1e-5
            assert orbit["max_sym_error_A_inv"] < 1e-5
            assert all(n["chirality"] == n["expected_chirality"] and n["isolation_eV"] > 1e-3
                       and abs(n["velocity_det"]) > 1e-5 for n in orbit["nodes"])
        assert max(covariance.values()) < 1e-5


if __name__ == "__main__":
    main()

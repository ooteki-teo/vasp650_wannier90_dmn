# Symmetry-adapted Wannier90 export for VASP 6.5.0

This extension adds the opt-in INCAR tag `LWRITE_DMN` to the VASP 6.5.0 Wannier90 interface. It exports the k-point stars, Bloch sewing matrices, and target-orbital representations required by Wannier90's `site_symmetry = .true.` workflow. The exporter supports scalar, collinear spin-polarized, and spinor calculations on a uniform k mesh closed under the selected unitary space-group operations. It includes PAW overlap contributions and the phases from nonsymmorphic translations and reciprocal-lattice folding.

The repository contains the extension, a recipe that generates a standard unified diff from a local licensed VASP 6.5.0 source tree, and regression inputs. The generated diff reproduces the previously tested source changes byte for byte. Start with [USAGE.md](USAGE.md) for installation, INCAR examples, and the validation commands.

## Validation

The extension was compiled from VASP 6.5.0 on macOS arm64 with GNU Fortran 14.2, Open MPI and OpenBLAS. Both CPU targets, `vasp_std` and `vasp_ncl`, were tested with standalone Wannier90 3.1.0.

| System | Result |
| --- | --- |
| Si | Four bond-centered Wannier functions. The 48 unitary operations produce k-point stars identical to the QE 7.5 reference, and Γ-point representation characters agree to the tested precision of 10⁻⁷. |
| AFM-II NiO | Both collinear channels complete Wannierization. The isolated O s case reaches the localization tolerance. The 16-orbital Ni d/O p case writes both Hamiltonians and reaches its disentanglement iteration limit. |
| SOC TaAs | Eight W1 and sixteen W2 nodes on the 8×8×8 model. Independently located nodes deviate from their predicted symmetry images by at most 6.64×10⁻⁶ Å⁻¹ for W1 and 6.04×10⁻⁷ Å⁻¹ for W2. The 24 chiralities follow the mirror and time-reversal relations. |

TaAs uses `symmetrize_eps = 1e-10`. The Weyl checks use the resulting Wannier Hamiltonian directly. These values measure symmetry in that interpolated model. Absolute node positions and energies depend on the electronic-structure and Wannierization settings. The [Weyl verification data](tests/reference/TaAs-weyl.json) include all 24 nodes.

## License

[MIT](LICENSE) covers the original material distributed in this repository. Applying the extension requires separately licensed VASP 6.5.0 and a compatible Wannier90 installation. VASP and the locally generated VASP sources remain subject to the applicable VASP license. [VASP installation documentation](https://vasp.at/wiki/Installing_VASP.6.X.X)

The file layout follows the [Wannier90 DMN specification](https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/files/#seednamedmn). The Si comparison uses [Quantum ESPRESSO's `write_dmn` interface](https://quantum-espresso.org/Doc/INPUT_pw2wannier90.html).

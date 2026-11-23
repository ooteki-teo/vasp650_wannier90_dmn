#!/usr/bin/env bash
set -euo pipefail

: "${VASP_STD:?Set VASP_STD to the absolute path to patched vasp_std}"
: "${VASP_NCL:?Set VASP_NCL to the absolute path to patched vasp_ncl}"
: "${WANNIER90:?Set WANNIER90 to the absolute path to wannier90.x}"
: "${POTCAR_ROOT:?Set POTCAR_ROOT to the VASP PAW-PBE library}"
: "${DMN_PYTHON:?Set DMN_PYTHON to Python with NumPy, SciPy, and spglib}"
if [[ $# != 1 ]]; then
    echo "Usage: $0 NEW_OUTPUT_DIRECTORY" >&2
    exit 2
fi

test_inputs="$(cd "$(dirname "$0")" && pwd)"
mkdir "$1"
cd "$1"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

copy_case() {
    local source="$1" dest="$2"
    mkdir "$dest"
    cp "$test_inputs/$source/"{INCAR,KPOINTS,POSCAR} "$dest/"
}

copy_case Si Si
copy_case NiO NiO
copy_case NiO-isolated NiO-isolated
copy_case TaAs-scf TaAs-scf
copy_case TaAs-final TaAs
cp "$POTCAR_ROOT/Si/POTCAR" Si/POTCAR
cat "$POTCAR_ROOT/Ni/POTCAR" "$POTCAR_ROOT/O/POTCAR" > NiO/POTCAR
ln -s ../NiO/POTCAR NiO-isolated/POTCAR
cat "$POTCAR_ROOT/Ta/POTCAR" "$POTCAR_ROOT/As/POTCAR" > TaAs-scf/POTCAR
for file in POTCAR CHGCAR WAVECAR; do
    ln -s "../TaAs-scf/$file" "TaAs/$file"
done

qe_check=(--check)
if [[ -n "${DMN_QE_BIN:-}" ]]; then
    : "${ESPRESSO_PSEUDO:?Set ESPRESSO_PSEUDO to the SSSP PBE directory for the optional QE check}"
    export ESPRESSO_PSEUDO
    mkdir Si-qe
    cp "$test_inputs/Si-qe/"{Si.win,scf.in,nscf.in,pw2wan.in} Si-qe/
    (
        cd Si-qe
        mpirun -np 2 "$DMN_QE_BIN/pw.x" -in scf.in > scf.log 2>&1
        mpirun -np 2 "$DMN_QE_BIN/pw.x" -in nscf.in > nscf.log 2>&1
        "$WANNIER90" -pp Si > w90-pp.log 2>&1
        mpirun -np 2 "$DMN_QE_BIN/pw2wannier90.x" -in pw2wan.in > pw2wan.log 2>&1
        "$WANNIER90" Si > w90.log 2>&1
        grep -q 'All done: wannier90 exiting' Si.wout
    )
    qe_check+=(--qe-seed Si-qe/Si)
fi

(
    cd Si
    mpirun -np 2 "$VASP_STD" > vasp.log 2>&1
    "$WANNIER90" wannier90 > w90.log 2>&1
    grep -q 'All done: wannier90 exiting' wannier90.wout
)
"$DMN_PYTHON" "$test_inputs/check_dmn.py" Si/wannier90 "${qe_check[@]}" > Si/check.json

run_nio() (
    local case_path="$1"
    shift
    cd "$case_path"
    mpirun -np 4 "$VASP_STD" > vasp.log 2>&1
    for spin in 1 2; do
        seed="wannier90.$spin"
        cp wannier90.win "$seed.win"
        "$WANNIER90" "$seed" > "w90.$spin.log" 2>&1
        grep -q 'All done: wannier90 exiting' "$seed.wout"
        "$DMN_PYTHON" "$test_inputs/check_dmn.py" "$seed" --magnetic --check "$@" > "check.$spin.json"
        if [[ "$case_path" == NiO-isolated ]]; then
            grep -q 'Wannierisation convergence criteria satisfied' "$seed.wout"
        fi
    done
)

run_nio NiO --window 20
run_nio NiO-isolated

(
    cd TaAs-scf
    mpirun -np 8 "$VASP_NCL" > vasp.log 2>&1
    test -s CHGCAR
    test -s WAVECAR
)
(
    cd TaAs
    mpirun -np 8 "$VASP_NCL" > vasp.log 2>&1
    "$WANNIER90" wannier90 > w90.log 2>&1
    grep -q 'All done: wannier90 exiting' wannier90.wout
    "$DMN_PYTHON" "$test_inputs/check_dmn.py" wannier90 --window 16.5 --check > check.json
    "$DMN_PYTHON" "$test_inputs/check_taas_weyl.py" wannier90 --check > weyl.json
)

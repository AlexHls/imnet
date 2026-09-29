#!/usr/bin/env python3
"""Compare a real imnet binary against the nine-species static-burn reference."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np

DATA = Path(__file__).resolve().parent / "data"


def prepare_network(source, destination, backend, species):
    """Copy mutable inputs; never run NuPPN in the supplied source directory."""
    destination.mkdir()
    if backend == "YANN":
        shutil.copyfile(DATA / "species9.txt", destination / "species.txt")
        for name in ("jinareaclib.dat", "part.txt", "mass.txt", "lmp_weak_rates.txt"):
            shutil.copyfile(source / name, destination / name)
        return

    npdata = (source.parent / "NPDATA").resolve(strict=True)
    (destination.parent / "NPDATA").symlink_to(npdata, target_is_directory=True)
    for name in ("ppn_frame.input", "ppn_solver.input", "ppn_physics.input",
                 "isotopedatabase.txt"):
        shutil.copyfile(source / name, destination / name)
    keep = {(a, z) for _, a, z in species}
    physics = destination / "ppn_physics.input"
    lines = physics.read_text().splitlines()
    found = set()
    for i, line in enumerate(lines):
        # NuPPN VITAL fixed format: I4,1X,A5,1X,2F4.0,1X,L1.
        if re.match(r"^\s*\d+ .{5} +\d+\. +\d+\. [TF]", line):
            az = (int(float(line[11:15])), int(float(line[15:19])))
            enabled = az in keep
            if enabled:
                found.add(az)
            lines[i] = line[:20] + ("T" if enabled else "F") + line[21:]
        elif re.match(r"^\s*\d+ [TF] ", line):
            # Use the backend's REACLIB rates, without VITAL overrides.
            lines[i] = line[:5] + "F" + line[6:]
    if found != keep:
        raise ValueError("NuPPN VITAL table does not contain all reference species")
    physics.write_text("\n".join(lines) + "\n")
    database = destination / "isotopedatabase.txt"
    lines = database.read_text().splitlines()
    for i in range(3, len(lines)):
        fields = lines[i].split()
        if len(fields) != 5:
            raise ValueError(f"Unexpected isotope database row: {lines[i]}")
        fields[4] = "T" if (int(fields[1]), int(fields[0])) in keep else "F"
        lines[i] = " ".join(fields)
    database.write_text("\n".join(lines) + "\n")


def run(command, cwd, log):
    with log.open("w") as output:
        result = subprocess.run([str(x) for x in command], cwd=cwd,
                                stdout=output, stderr=subprocess.STDOUT, timeout=240)
    if result.returncode:
        raise RuntimeError(f"Run failed ({result.returncode}); see {log}")


def load_state(path, backend, species, times):
    # Some Fortran STOP paths return zero: require a complete, fresh export too.
    state = json.loads(path.read_text())
    assert state["network_backend"].upper() == backend
    assert state["units"]["abundance"] == "mass_fraction"
    metadata = state["species"]
    names = [s["name"] for s in metadata]
    expected = [s[0] for s in species]
    assert len(names) == len(expected) and set(names) == set(expected), names
    order = [names.index(name) for name in expected]
    for name, a, z in species:
        s = metadata[names.index(name)]
        assert (s["A"], s["Z"], s["N"]) == (a, z, a - z)
    steps = state["steps"]
    np.testing.assert_allclose([s["time"] for s in steps], times, rtol=1e-13, atol=0)
    assert all(s["success"] and s["rho"] == 1e4 and s["temp"] == 2e8 for s in steps)
    values = np.array([s["abundances"] for s in steps])[:, order]
    assert np.isfinite(values).all() and (values >= 0).all()
    np.testing.assert_allclose(values.sum(axis=1), 1, rtol=0, atol=1e-6)
    np.testing.assert_allclose(values[0], [0.5, 0.25, 0.25, 0, 0, 0, 0, 0, 0], atol=1e-15)
    edges = set()
    for step in steps:
        assert np.isfinite(step["dedt"])
        for flux in step["fluxes"]:
            src, dst = flux["source"], flux["target"]
            assert 0 <= src < len(names) and 0 <= dst < len(names) and src != dst
            for metric in ("strength_rate", "strength_dydt", "strength_dxdt"):
                assert np.isfinite(flux[metric]) and flux[metric] > 0
            assert flux["equation"]
            edges.add((names[src], names[dst], flux["weak"]))
    # Expected directed hot-CNO flow, including beta-decays.
    required = {("c12", "n13", False), ("n13", "c13", True),
                ("c13", "n14", False), ("n13", "o14", False),
                ("o14", "n14", True), ("n14", "o15", False),
                ("o15", "n15", True), ("n15", "c12", False)}
    assert required <= edges, f"Missing flux arrows: {required - edges}"
    return state, values


def plot_results(output, times, reference, results, species, rtol, atol, state):
    os.environ.setdefault("MPLCONFIGDIR", str(output / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(output / ".cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 3, figsize=(13, 10), sharex=True)
    for i, (ax, (name, _, _)) in enumerate(zip(axes.flat, species)):
        ax.loglog(times[1:], np.maximum(reference[1:, i], 1e-12), "k-", label="Reference")
        for label, values in results.items():
            ax.loglog(times[1:], np.maximum(values[1:, i], 1e-12), "--", label=label)
        ax.set(title=name, ylim=(1e-12, 1), ylabel="Mass fraction X")
        ax.grid(True, alpha=0.2)
    axes[0, 0].legend(fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)")
    fig.suptitle("Static burn: T = 2e8 K, rho = 1e4 g/cm³; display floor 1e-12")
    fig.tight_layout()
    fig.savefig(output / "abundances.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    for label, values in results.items():
        error = np.max(abs(values - reference) / (atol + rtol * reference), axis=1)
        ax.loglog(times[1:], np.maximum(error[1:], 1e-12), label=label)
    ax.axhline(1, color="black", linestyle=":", label="Pass limit")
    ax.set(xlabel="Time (s)", ylabel="Worst species error / tolerance",
           title=f"Every sample must satisfy |X − reference| ≤ {atol:g} + {rtol:g} × reference")
    ax.legend()
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(output / "errors.png", dpi=160)
    plt.close(fig)


    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharex=True, sharey=True)
    metadata = state["species"]
    for ax, target_time in zip(axes, (1, 100, 1000)):
        step = min(state["steps"], key=lambda s: abs(s["time"] - target_time))
        # Heavy-nucleus CNO arrows only; omit p/alpha legs to avoid clutter.
        fluxes = [f for f in step["fluxes"] if metadata[f["source"]]["A"] >= 12
                  and metadata[f["target"]]["A"] >= 12 and f["strength_dydt"] > 1e-12]
        maximum = max((f["strength_dydt"] for f in fluxes), default=1)
        for f in fluxes:
            a, b = metadata[f["source"]], metadata[f["target"]]
            ax.annotate("", xy=(b["N"], b["Z"]), xytext=(a["N"], a["Z"]),
                        arrowprops=dict(arrowstyle="->", color="tab:orange" if f["weak"] else "tab:blue",
                                        lw=0.5 + 3 * np.sqrt(f["strength_dydt"] / maximum),
                                        shrinkA=15, shrinkB=15, connectionstyle="arc3,rad=0.08"))
        for sp in metadata:
            if sp["A"] >= 12:
                ax.text(sp["N"], sp["Z"], sp["name"], ha="center", va="center",
                        bbox=dict(boxstyle="round", fc="white", ec="0.7"))
        ax.set(title=f"t = {step['time']:.3g} s", xlabel="Neutron number N",
               xlim=(4.5, 8.5), ylim=(5.5, 8.5), xticks=range(5, 9), yticks=range(6, 9))
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Proton number Z")
    fig.suptitle("Exported CNO fluxes: blue = regular, orange = weak; width ∝ √|dY/dt| (per panel)")
    fig.tight_layout()
    fig.savefig(output / "fluxes.png", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--imnet", type=Path, required=True)
    parser.add_argument("--backend", choices=["NUPPN", "YANN"], required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--gui-driver", type=Path, help="AppState cached-step/flux test executable")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rtol", type=float, default=0.05)
    parser.add_argument("--atol", type=float, default=1e-5)
    args = parser.parse_args()
    if not np.isfinite([args.rtol, args.atol]).all() or min(args.rtol, args.atol) <= 0:
        parser.error("Tolerances must be finite and positive")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    species = [(name, int(a), int(z)) for name, a, z in
               (line.split() for line in (DATA / "species9.txt").read_text().splitlines()[1:] if line.strip())]
    with np.load(DATA / "X_reference.npz", allow_pickle=False) as data:
        times, reference = data["t"], data["X"]
    assert reference.shape == (len(times), len(species))
    assert times[0] == 0 and np.isfinite(times).all() and (np.diff(times) > 0).all()
    # The supplied solver reference has roundoff negatives down to -1.5e-17.
    assert np.isfinite(reference).all() and (reference >= -1e-15).all()
    reference = np.maximum(reference, 0)
    np.testing.assert_allclose(reference.sum(axis=1), 1, atol=1e-6, rtol=0)
    # Preserve each isolated run (including its inputs and logs) for GUI inspection.
    work = Path(tempfile.mkdtemp(prefix="run-", dir=output))
    network = work / "network"
    prepare_network(args.data_dir.resolve(), network, args.backend, species)
    initial, trajectory = work / "initial.txt", work / "trajectory.txt"
    initial.write_text("p 0.5\nhe4 0.25\nc12 0.25\n")
    np.savetxt(trajectory, np.column_stack([times, np.full_like(times, 1e4), np.full_like(times, 2e8)]), fmt="%.17e")
    # CLI failures must preserve inputs and produce a nonzero status.
    binary = args.imnet.resolve()
    for extra in (["--output", initial], ["--output", trajectory],
                  ["--output", work / "same", "--save-state", work / "same"],
                  ["--output", "-"]):
        completed = subprocess.run([str(x) for x in
            [binary, "--headless", "--data-dir", network, "--abundances", initial,
             "--trajectory", trajectory, *extra]], cwd=work, capture_output=True, timeout=20)
        assert completed.returncode != 0, "Invalid CLI invocation was accepted"
    assert initial.read_text() == "p 0.5\nhe4 0.25\nc12 0.25\n"
    assert len(np.loadtxt(trajectory)) == len(times)
    single = work / "single.json"
    run([binary, "--headless", "--data-dir", network, "--abundances", initial,
         "--rho", "1e4", "--temp", "2e8", "--dt", "1e-6", "--save-state", single],
        work, work / "single.log")
    single_state = json.loads(single.read_text())
    assert single_state["steps"][0]["time"] == 1e-6
    zero = work / "zero.txt"
    zero.write_text("p 0\n")
    completed = subprocess.run([str(x) for x in
        [binary, "--headless", "--data-dir", network, "--abundances", zero,
         "--output", work / "zero.csv"]], cwd=work, capture_output=True, timeout=30)
    assert completed.returncode != 0, "All-zero composition was accepted"
    state_path, csv_path = work / "headless.json", work / "headless.csv"
    run([args.imnet.resolve(), "--headless", "--data-dir", network, "--abundances", initial,
         "--trajectory", trajectory, "--output", csv_path, "--save-state", state_path],
        work, work / "headless.log")
    state, values = load_state(state_path, args.backend, species, times)
    with csv_path.open() as file:
        rows = list(csv.DictReader(file))
    np.testing.assert_allclose([float(row["time"]) for row in rows], times, rtol=1e-13, atol=0)
    np.testing.assert_allclose([[float(row[name]) for name, _, _ in species] for row in rows],
                               values, rtol=1e-13, atol=0)
    results = {f"{args.backend} headless": values}
    if args.gui_driver:
        gui_path = work / "gui-state.json"
        run([args.gui_driver.resolve(), network, initial, trajectory, gui_path], work, work / "gui-state.log")
        _, gui = load_state(gui_path, args.backend, species, times)
        results[f"{args.backend} GUI state"] = gui
    plot_results(output, times, reference, results, species, args.rtol, args.atol, state)
    summary = {"backend": args.backend, "rtol": args.rtol, "atol": args.atol,
               "inputs_and_logs": str(work), "results": {}}
    passed = True
    for label, actual in results.items():
        scaled = abs(actual - reference) / (args.atol + args.rtol * reference)
        i, j = np.unravel_index(scaled.argmax(), scaled.shape)
        ok = bool((scaled <= 1).all())
        passed &= ok
        summary["results"][label] = {"passed": ok, "worst_scaled_error": float(scaled[i, j]),
            "time": float(times[i]), "species": species[j][0],
            "reference": float(reference[i, j]), "actual": float(actual[i, j]),
            "max_absolute_error_by_species": dict(zip([s[0] for s in species],
                                                       abs(actual-reference).max(axis=0).tolist()))}
    if args.gui_driver:
        np.testing.assert_allclose(gui, values, rtol=1e-12, atol=1e-15,
                                   err_msg="GUI state and headless disagree")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Plots: {output / 'abundances.png'}, {output / 'errors.png'}, {output / 'fluxes.png'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

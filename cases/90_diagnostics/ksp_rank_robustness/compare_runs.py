"""Compare the final flat frame of two pancake run directories, field by field."""
import sys
from pathlib import Path

sys.path.insert(0, "/tmp/ksp_repro")
from solver_io import final_pvts, _read_pvts  # noqa: E402


def load(run_dir):
    p = final_pvts(run_dir)
    if p is None:
        raise SystemExit(f"no frame in {run_dir}")
    return {k: v[0] for k, v in _read_pvts(p).items()}


def main(a, b):
    fa, fb = load(a), load(b)
    print(f"{'field':<24}{'n':>8}{'max|d|':>14}{'max|d|/range':>16}")
    worst = 0.0
    for name in sorted(set(fa) & set(fb)):
        va, vb = fa[name], fb[name]
        if len(va) != len(vb):
            print(f"{name:<24}{'LEN MISMATCH':>8}")
            continue
        d = max(abs(x - y) for x, y in zip(va, vb))
        rng = max(va) - min(va)
        rel = d / rng if rng > 0 else 0.0
        worst = max(worst, rel)
        print(f"{name:<24}{len(va):>8}{d:>14.4e}{rel:>16.4e}")
    print(f"worst relative deviation: {worst:.4e}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))

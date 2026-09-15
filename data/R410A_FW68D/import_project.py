"""Decode chart pixels to SI-ready CSV, using only the Python standard library.

Usage: python import_project.py [fw68d-project.json]
Without an argument, regenerate from the archived project without its image.
Pressure is absolute and concentration is refrigerant mass / solution mass,
as specified by the user, overriding the project's unspecified basis labels.
"""
import bisect
import csv
import json
import math
from pathlib import Path
import sys


def decode_axis(pixel, axis):
    """Interpolate adjacent calibration anchors in the declared transformed axis."""
    anchors = sorted(axis["anchors"], key=lambda item: item["pixel"])
    pixels = [float(item["pixel"]) for item in anchors]
    if len(pixels) < 2 or any(a >= b for a, b in zip(pixels, pixels[1:])):
        raise ValueError("Calibration needs distinct pixel anchors")
    scale = axis["scale"]
    if scale not in ("linear", "log", "logK"):
        raise ValueError(f"Unsupported calibration scale: {scale}")
    offset = 273.15 if scale == "logK" else 0.0
    values = [float(item["value"]) + offset for item in anchors]
    if scale != "linear":
        values = [math.log(value) for value in values]
    i = min(max(bisect.bisect_right(pixels, pixel), 1), len(pixels) - 1)
    weight = (pixel - pixels[i - 1]) / (pixels[i] - pixels[i - 1])
    value = values[i - 1] + weight * (values[i] - values[i - 1])
    return (value if scale == "linear" else math.exp(value)) - offset


def convert(project):
    """Return every concentration/oil point; pressure-free oil rows have blank P."""
    rows = []
    for curve in project["curves"]:
        if curve["kind"] not in ("concentration", "oil"):
            continue
        concentration = float(curve["value"]) / 100.0
        for point in curve["points"]:
            temperature = decode_axis(float(point["x"]), project["calibration"]["x"])
            viscosity = decode_axis(float(point["y"]), project["calibration"]["y"])
            pressure = "" if curve["kind"] == "oil" else float(point["p_MPa"])
            if not (0 <= concentration <= 1 and temperature > -273.15 and viscosity > 0):
                raise ValueError("Nonphysical chart point")
            if pressure != "" and not (math.isfinite(pressure) and pressure > 0):
                raise ValueError("Pressure must be positive absolute MPa")
            if not all(math.isfinite(v) for v in (concentration, temperature, viscosity)):
                raise ValueError("Nonfinite chart point")
            rows.append((concentration, temperature, pressure, viscosity))
    return sorted(rows)


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else directory / "source_project.json"
    project = json.loads(source.read_text(encoding="utf-8-sig"))
    rows = convert(project)
    with (directory / "fw68d.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("c_mass_fraction", "T_C", "P_MPa", "viscosity_mm2s"))
        writer.writerows(rows)
    print(f"Decoded {len(rows)} points to {directory / 'fw68d.csv'}")

"""Run with python -m unittest discover -s data/R410A_FW68D."""
import csv
import json
from pathlib import Path
import unittest

from import_project import convert, decode_axis


class ImportProjectTest(unittest.TestCase):
    def test_calibration_and_complete_export(self):
        directory = Path(__file__).resolve().parent
        project = json.loads((directory / "source_project.json").read_text(encoding="utf-8-sig"))
        for axis in project["calibration"].values():
            if not isinstance(axis, dict):
                continue
            for anchor in axis["anchors"]:
                self.assertAlmostEqual(decode_axis(anchor["pixel"], axis), anchor["value"], places=7)
        self.assertAlmostEqual(decode_axis(0.5, {"scale": "log", "anchors": [
            {"pixel": 0, "value": 1}, {"pixel": 1, "value": 100}]}), 10)
        rows = convert(project)
        self.assertEqual(len(rows), 125)
        self.assertEqual(sum(row[0] == 0 for row in rows), 15)
        self.assertEqual(sum(row[0] == 1 for row in rows), 11)
        with (directory / "fw68d.csv").open(newline="", encoding="utf-8") as stream:
            saved = list(csv.reader(stream))[1:]
        self.assertEqual([[str(value) for value in row] for row in rows], saved)


if __name__ == "__main__":
    unittest.main()

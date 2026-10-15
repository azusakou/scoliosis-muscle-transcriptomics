"""Output guards use small valid blank-PDF fixtures, not scientific figures."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("figure_guard", ROOT / "scripts/verify_replayed_figures.py")
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


def blank_pdf(width, height):
    data = b"%PDF-1.4\n" + b"% parser fixture padding\n" * 60
    objects = [b"<< /Type /Catalog /Pages 2 0 R >>", b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
               (f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width*72/25.4} {height*72/25.4}] /Contents 4 0 R >>").encode(),
               b"<< /Length 0 >>\nstream\n\nendstream"]
    offsets = [0]
    for number, value in enumerate(objects, 1):
        offsets.append(len(data))
        data += f"{number} 0 obj\n".encode() + value + b"\nendobj\n"
    start = len(data)
    data += b"xref\n0 5\n0000000000 65535 f \n"
    for offset in offsets[1:]:
        data += f"{offset:010d} 00000 n \n".encode()
    return data + f"trailer\n<< /Size 5 /Root 1 0 R >>\nstartxref\n{start}\n%%EOF\n".encode()


class ReplayedFigureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.config = json.loads((ROOT / guard.CONFIG).read_text())
        path = self.root / guard.CONFIG
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(self.config))
        for item in self.config["figures"]:
            p = self.root / item["producer"]
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("# fixture producer identity\n")
        self.begun = guard.start(self.root)
        for item in self.config["figures"]:
            p = self.root / item["path"]
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(blank_pdf(item["width_mm"], item["height_mm"]))
            os.utime(p, ns=(self.begun["started_ns"]+1000000, self.begun["started_ns"]+1000000))

    def tearDown(self):
        self.temp.cleanup()

    def test_complete_fresh_outputs_have_hashes(self):
        result = guard.verify(self.root)
        self.assertEqual(result["generated_PDFs"], 16)
        self.assertTrue(all(len(r["sha256"]) == 64 for r in result["outputs"]))

    def test_missing_or_stale_output_is_rejected(self):
        p = self.root / self.config["figures"][0]["path"]
        os.utime(p, ns=(1, 1))
        with self.assertRaisesRegex(ValueError, "Stale/copied"):
            guard.verify(self.root)
        p.unlink()
        with self.assertRaisesRegex(ValueError, "Missing or empty"):
            guard.verify(self.root)

    def test_non_pdf_or_bad_xref_is_rejected(self):
        p = self.root / self.config["figures"][0]["path"]
        data = p.read_bytes().replace(b"startxref\n", b"startxref\n999999")
        p.write_bytes(data)
        os.utime(p, ns=(self.begun["started_ns"]+1000000, self.begun["started_ns"]+1000000))
        with self.assertRaisesRegex(ValueError, "Invalid/incomplete"):
            guard.verify(self.root)

    def test_supplied_pdf_cannot_be_counted_as_generated(self):
        (self.root / "SHA256SUMS.tsv").write_text("path\tbytes\tsha256\n" + self.config["figures"][0]["path"] + "\t1\tfixture\n")
        with self.assertRaisesRegex(ValueError, "immutable release inputs"):
            guard.verify(self.root)

    def test_changed_producer_is_rejected(self):
        (self.root / self.config["figures"][0]["producer"]).write_text("# changed producer\n")
        with self.assertRaisesRegex(ValueError, "producer changed"):
            guard.verify(self.root)


if __name__ == "__main__":
    unittest.main()

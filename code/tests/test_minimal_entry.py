"""The single entry must preserve strings and reject unsafe/corrupt inputs."""
import csv,gzip,importlib.util,json,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('minimal_entry',ROOT/'run.py')
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)

class MinimalEntryTests(unittest.TestCase):
    def fixture(self,d):
        root=Path(d);p=root/'data/input.tsv';p.parent.mkdir();p.write_text('gene\tcount\nG\t2\n')
        with (root/'MANIFEST.tsv').open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['path','runtime','bytes','sha256']);w.writerow(['data/input.tsv','input.tsv',p.stat().st_size,r.reproduction.digest(p)])
        return root,p
    def test_manifest_hash_and_complete_inventory(self):
        with tempfile.TemporaryDirectory() as d:
            root,p=self.fixture(d);self.assertEqual(len(r.verify(root)),1)
            p.write_text('gene\tcount\nG\t3\n')
            with self.assertRaises(ValueError):r.verify(root)
        with tempfile.TemporaryDirectory() as d:
            root,p=self.fixture(d);(root/'extra').write_text('untracked')
            with self.assertRaisesRegex(ValueError,'Uninventoried'):r.verify(root)
    def test_table_roundtrip_preserves_exact_text(self):
        with tempfile.TemporaryDirectory() as d:
            table={'columns':['gene','number','text'],'rows':[['001','0.050000000000000044','tab\tquote"\nline']]}
            r.restore_table(Path(d),'analysis/a.tsv.gz',table)
            with gzip.open(Path(d)/'analysis/a.tsv.gz','rt',newline='') as f:rows=list(csv.reader(f,delimiter='\t'))
            self.assertEqual(rows,[table['columns'],*table['rows']])
            with self.assertRaises(ValueError):r.restore_table(Path(d),'analysis/a.tsv.gz',table)
    def test_unsafe_paths_bad_structure_and_numeric_json_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for name in ['../escape.tsv','/tmp/a.tsv','analysis//a.tsv','analysis\\a.tsv']:
                with self.subTest(path=name),self.assertRaises(ValueError):r.restore_table(root,name,{'columns':['a'],'rows':[['1']]})
            for table in [{'columns':['a'],'rows':[[1]]},{'columns':['a'],'rows':[['1','2']]},{'columns':['a','a'],'rows':[]}]:
                with self.assertRaises(ValueError):r.restore_table(root,'analysis/b.tsv',table)
    def test_symlink_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root,p=self.fixture(d);target=root/'target';p.rename(target);p.symlink_to(target)
            with self.assertRaisesRegex(ValueError,'Symlink'):r.verify(root)
if __name__=='__main__':unittest.main()

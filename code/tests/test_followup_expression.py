import importlib,sys,unittest
from decimal import Decimal
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
m=importlib.import_module('24_export_followup_gene_expression')

class FollowupExpressionTests(unittest.TestCase):
    def test_gene_sums_preserve_exact_decimal_ties_and_mapping(self):
        rows=[dict(Name='a',TPM='0.100000',NumReads='1.001'),
              dict(Name='b',TPM='0.200000',NumReads='2.002'),
              dict(Name='c',TPM='0.300000',NumReads='3.003'),
              dict(Name='unmapped',TPM='0.000000',NumReads='0.000')]
        result,mapped,total=m.aggregate_gene_rows(rows,{'a':'gene1','b':'gene1','c':'gene2'})
        self.assertEqual(result['gene1'],result['gene2'])
        self.assertEqual(result['gene1'],[Decimal('.300000'),Decimal('3.003')])
        self.assertEqual((mapped,total),(3,4))
    def test_nonfinite_negative_and_malformed_numbers_are_rejected(self):
        for value in ['NaN','Infinity','-0.1','invalid']:
            with self.assertRaises(ValueError):m.aggregate_gene_rows([dict(Name='a',TPM=value,NumReads='1')],{'a':'g'})

if __name__=='__main__':unittest.main()

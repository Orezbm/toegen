import unittest
from toegen.homology import HomologySwitchGenerator

class TestHomologySwitchGenerator(unittest.TestCase):

    def setUp(self):
        # Example genome data to simulate transcript sequences
        self.transcripts_dict = [
            {'gene': 'GeneA', 'protein': 'ProteinA', 'sequence': 'AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG'},
            {'gene': 'GeneB', 'protein': 'ProteinB', 'sequence': 'AUGGCUCGAUUGCCCUUCUAGGAUCGUAGCUAGGAUC'}
        ]
        self.reporter_gene = "AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG"
        self.homology_generator = HomologySwitchGenerator(
            cell_type="Homo sapiens",
            reporter_gene=self.reporter_gene,
            transcripts_dict=self.transcripts_dict
        )

    def test_generate_switch_with_homology(self):
        triggers = ["AUGGCCUAGCGCUAUGCCCUAGG", "AUGGCUCGAUUGCCCUUCUAGGAUCGUAGC"]  # Eukaryotic and Prokaryotic
        result_df = self.homology_generator.generate_switch_with_homology(triggers, self.transcripts_dict)
        self.assertEqual(len(result_df), 2)  # Two triggers processed
        self.assertTrue("switch" in result_df.columns)
        self.assertTrue("complex_concentration" in result_df.columns)

    def test_homology_search(self):
        triggers = ["AUGGCCUAGCGCUAUGCCCUAGG"]
        homology_results = self.homology_generator.find_homology(triggers[0], self.transcripts_dict)
        self.assertGreaterEqual(len(homology_results), 1)

    def test_invalid_trigger(self):
        # Invalid length trigger
        triggers = ["AUGGCUCGAU"]  # Short trigger, less than 23 nucleotides
        with self.assertRaises(ValueError):
            self.homology_generator.generate_switch_with_homology(triggers, self.transcripts_dict)


if __name__ == "__main__":
    unittest.main()

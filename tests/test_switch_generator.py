import unittest
from toegen.switch_generator import SwitchGenerator, ProkaryoticSwitchGenerator


class TestSwitchGenerator(unittest.TestCase):

    def setUp(self):
        self.translated_gene = "AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG"
        self.euk_switch_generator = SwitchGenerator(translated_gene_sequence=self.translated_gene)
        self.prok_switch_generator = ProkaryoticSwitchGenerator(translated_gene_sequence=self.translated_gene)

    def test_eukaryotic_trigger_length(self):
        trigger = "AUGGCCUAGCGCUAUGCCCUAGG"  # 23 nucleotides for eukaryotic
        switch_strand, defect, concentration, mfe_structure = self.euk_switch_generator.get_switch(trigger)
        self.assertEqual(len(trigger), 23)
        self.assertTrue(isinstance(switch_strand, str))
        self.assertTrue(isinstance(defect, float))
        self.assertTrue(isinstance(concentration, float))

    def test_prokaryotic_trigger_length(self):
        trigger = "AUGGCUCGAUUGCCCUUCUAGGAUCGUAGC"  # 30 nucleotides for prokaryotic
        switch_strand, defect, concentration, mfe_structure = self.prok_switch_generator.get_switch(trigger)
        self.assertEqual(len(trigger), 30)
        self.assertTrue(isinstance(switch_strand, str))
        self.assertTrue(isinstance(defect, float))
        self.assertTrue(isinstance(concentration, float))

    def test_invalid_trigger(self):
        # Invalid trigger length
        with self.assertRaises(ValueError):
            self.euk_switch_generator.get_switch("AUGGCCUAG")  # Only 9 nucleotides, should raise error

        # Invalid characters
        with self.assertRaises(ValueError):
            self.euk_switch_generator.get_switch("AUGGXXUAGCGCUAUGCCCUAGG")  # Contains 'X'


if __name__ == "__main__":
    unittest.main()

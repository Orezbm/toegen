from Bio.Seq import Seq
from nupack import *
from typing import Tuple


config.threads = 8

KOZAK_SEQUENCE = "ACCAAAAUG"
RBS_SEQUENCE = "AGAGGAGA"
STOP_CODONS = ['UGA', 'UAG', 'UAA']
FORBIDDEN_PATTERNS = ["AAAA", "CCCC", "GGGG", "UUUU", "KKKKKK", "MMMMMM", "RRRRRR", "SSSSSS", "WWWWWW", "YYYYYY"]
START_CODON = "AUG"

class SwitchGenerator:

    SEED = 42
    VALID_NUCLEOTIDES = {'A', 'C', 'G', 'U'}

    def __init__(self, translated_gene_sequence: str, leader_sequence: str = "",
                 a_domain_size: int = 8, b_domain_size: int = 15, extra_stem_length: int = 0,
                 loop_size: int = 12, linker_pattern: str = "",
                 nucleotides_count_to_trim_gene_after: int = 50, is_kozak_down: bool = True):
        self.model = Model(material='rna', celsius=37)
        self.is_kozak_down = is_kozak_down
        self.translated_gene_sequence = translated_gene_sequence[3:nucleotides_count_to_trim_gene_after + 3]
        self.b_domain_size = b_domain_size
        self.a_domain_size = a_domain_size
        self.extra_stem_length = extra_stem_length
        self.loop_size = loop_size
        self.trigger_length = a_domain_size + b_domain_size
        self.leader = Domain(leader_sequence, name="leader", material="rna")
        self.extra_stem = Domain(f"N{extra_stem_length}", name="b end", material="rna")
        self.linker = Domain(linker_pattern, name="linker", material="rna")
        self.translated_gene = Domain(self.translated_gene_sequence, name="translated gene", material="rna")

        if is_kozak_down:
            self.loop = Domain(f"N{loop_size}", name="loop", material="rna")
            self.kozak = Domain(f"N{KOZAK_SEQUENCE}", name="kozak", material="rna")
            self.switch_structure = (f".{len(leader_sequence)}.{a_domain_size}({b_domain_size + extra_stem_length}"
                                     f".{loop_size}){b_domain_size + extra_stem_length}.{len(KOZAK_SEQUENCE) + 1}"
                                     f".{len(linker_pattern)}.{len(self.translated_gene_sequence)}")
        else:
            assert loop_size > len(KOZAK_SEQUENCE), "In Kozak up mode - loop size must be larger than the kozak length"
            self.loop = Domain(f"N{loop_size - len(KOZAK_SEQUENCE)}{KOZAK_SEQUENCE}", name="loop", material="rna")
            post_stem_offset = 3 - ((b_domain_size + extra_stem_length + len(linker_pattern)) % 3)
            self.post_stem = Domain(f"N{post_stem_offset}", name="post stem", material="rna")
            self.switch_structure = (f".{len(leader_sequence)}.{a_domain_size}({b_domain_size + extra_stem_length}"
                                     f".{loop_size}){b_domain_size + extra_stem_length}.{post_stem_offset}"
                                     f".{len(linker_pattern)}.{len(self.translated_gene_sequence)}")

        self.patterns_to_prevent = Pattern(FORBIDDEN_PATTERNS)
        self.options = DesignOptions(f_stop=0.02, seed=self.SEED, wobble_mutations=True)
        self.complex_name = 'activated_complex'

    @property
    def prescribed_stem_structure(self):
        return (f"({self.b_domain_size + self.extra_stem_length}"
                f".{self.loop_size}){self.b_domain_size + self.extra_stem_length}")

    def get_switch(self, trigger_sequence: str, healthy_sequence: str = None, is_watson_crick: bool = False) -> Tuple[
        str, float, float, str]:

        if len(trigger_sequence) != self.trigger_length:
            raise ValueError(f"The trigger sequence must be exactly {self.trigger_length} nucleotides long for the eukaryotic/human switch, but you provided {len(trigger_sequence)} nucleotides.")

        if not set(trigger_sequence).issubset(self.VALID_NUCLEOTIDES):
            invalid_chars = set(trigger_sequence) - self.VALID_NUCLEOTIDES
            raise ValueError(f"The trigger sequence contains invalid characters: {invalid_chars}. Only 'A', 'C', 'G', and 'U' are allowed.")

        assert len(trigger_sequence) == self.trigger_length

        trigger_domain = Domain(trigger_sequence, name="trigger domain", material='rna')
        trigger = TargetStrand([trigger_domain], name="trigger")
        trigger_structure = "." * len(trigger_sequence)
        trigger_complex = TargetComplex([trigger], trigger_structure, name='trigger complex')

        if is_watson_crick:
            trigger_reverse_complement = str(Seq(trigger_sequence).reverse_complement_rna())
            a_domain = Domain(f"{trigger_reverse_complement[:self.a_domain_size]}", name="a", material="rna")
            b_domain = Domain(f"{trigger_reverse_complement[-self.b_domain_size:]}", name="b", material="rna")
        else:
            a_domain = Domain(f"N{self.a_domain_size}", name="a", material="rna")
            b_domain = Domain(f"N{self.b_domain_size}", name="b", material="rna")

        if self.is_kozak_down:
            switch = TargetStrand([
                self.leader, a_domain, b_domain, self.extra_stem,
                self.loop, ~self.extra_stem, ~b_domain, self.kozak, self.linker, self.translated_gene
            ], name="switch")
        else:
            switch = TargetStrand([
                self.leader, a_domain, b_domain, self.extra_stem,
                self.loop, ~self.extra_stem, ~b_domain, self.post_stem, self.linker, self.translated_gene
            ], name="switch")
            stop_codon_pattern = Pattern(STOP_CODONS, scope=[~self.extra_stem, ~b_domain, self.post_stem])

        switch_complex = TargetComplex([switch], self.switch_structure, name='switch complex')

        if healthy_sequence:
            healthy_sequence_domain = Domain(healthy_sequence, name="healthy trigger domain", material='rna')
            healthy_trigger = TargetStrand([healthy_sequence_domain], name="healthy trigger")
            healthy_trigger_structure = "." * len(healthy_sequence)
            healthy_trigger_complex = TargetComplex([healthy_trigger], healthy_trigger_structure,
                                                    name='healthy trigger complex')

        activated_structure = self._get_activated_structure(switch)
        activated_complex = TargetComplex([trigger, switch], activated_structure, name=self.complex_name)

        reactants = TargetTube(on_targets={switch_complex: 1e-08, trigger_complex: 1e-08},
                               off_targets=SetSpec(max_size=2, exclude=tuple([activated_complex])), name='reactants')

        if healthy_sequence:
            off_target_reactants = TargetTube(on_targets={switch_complex: 1e-08, healthy_trigger_complex: 1e-08},
                                              off_targets=SetSpec(max_size=2), name='off_target_reactants')

        products = TargetTube(on_targets={activated_complex: 1e-08}, off_targets=SetSpec(max_size=2), name='products')

        if healthy_sequence:
            tubes = [reactants, products, off_target_reactants]
        else:
            tubes = [reactants, products]

        my_design = tube_design(tubes=tubes,
                                hard_constraints=[] if is_watson_crick or self.is_kozak_down else [stop_codon_pattern],
                                soft_constraints=[],
                                model=self.model, options=self.options)

        design_results = my_design.run(trials=1)[0]

        switch_designed_strand = [v for k, v in design_results.to_analysis.strands.items() if
                                  v.name == "switch"].pop()
        ensemble_defect = self._get_ensemble_defect(design_results)
        complex_concentration = self._get_trigger_switch_complex_concentration(design_results, self.complex_name)
        switch_mfe_structure = self._get_switch_mfe_structure(design_results, switch_complex)
        return str(switch_designed_strand), ensemble_defect, complex_concentration, switch_mfe_structure

    def _get_activated_structure(self, switch) -> str:
        return ('(' * self.trigger_length + '+' + '.' * self.leader.nt() + ')' * self.trigger_length +
                '.' * (switch.nt() - self.leader.nt() - self.trigger_length))

    def _get_ensemble_defect(self, design_result) -> float:
        return design_result.defects.ensemble_defect

    def _get_trigger_switch_complex_concentration(self, design_result, complex_name) -> float:
        complex_concentrations = design_result.concentrations.table[
            design_result.concentrations.table.complex_name == complex_name].iloc[0]
        target_concentration = complex_concentrations["target_concentration"]
        predicted_concentration = complex_concentrations["concentration"]
        return predicted_concentration / target_concentration

    def _get_switch_mfe_structure(self, design_result, switch_target_strand) -> str:
        return str(mfe(strands=design_result.to_analysis(switch_target_strand), model=self.model)[0].structure)


class ProkaryoticSwitchGenerator:

    SEED = 42
    VALID_NUCLEOTIDES = {'A', 'C', 'G', 'U'}

    def __init__(self, translated_gene_sequence: str, leader_sequence: str = "GGG",
                 a_domain_size: int = 12, b_domain_size_pre_bulge: int = 9, b_domain_size_post_bulge: int = 6,
                 loop_size: int = 11, linker_pattern: str = "AACCUGGCGGCAGCGCAAAAG",
                 nucleotides_count_to_trim_gene_after: int = 30):
        self.model = Model(material='rna', celsius=37)
        self.translated_gene_sequence = translated_gene_sequence[3:nucleotides_count_to_trim_gene_after + 3]
        self.b_domain_size_pre_bulge = b_domain_size_pre_bulge
        self.b_domain_size_post_bulge = b_domain_size_post_bulge
        self.a_domain_size = a_domain_size
        self.loop_size = loop_size
        self.trigger_length = a_domain_size + b_domain_size_pre_bulge + b_domain_size_post_bulge + len(START_CODON)
        self.leader = Domain(leader_sequence, name="leader", material="rna")
        self.linker = Domain(linker_pattern, name="linker", material="rna")
        self.translated_gene = Domain(self.translated_gene_sequence, name="translated gene", material="rna")

        assert loop_size > len(RBS_SEQUENCE), "loop size must be larger than the RBS length"

        self.loop = Domain(f"N{loop_size - len(RBS_SEQUENCE)}{RBS_SEQUENCE}", name="loop", material="rna")
        self.bulge = Domain(f"{START_CODON}", name="bulge", material="rna")
        self.switch_structure = (f".{len(leader_sequence)}.{a_domain_size}({b_domain_size_pre_bulge}"
                                 f".{len(START_CODON)}({b_domain_size_post_bulge}"
                                 f".{loop_size}){b_domain_size_post_bulge}.{len(START_CODON)}"
                                 f"){b_domain_size_pre_bulge}"
                                 f".{len(linker_pattern)}.{len(self.translated_gene_sequence)}")

        self.patterns_to_prevent = Pattern(FORBIDDEN_PATTERNS)
        self.options = DesignOptions(f_stop=0.02, seed=self.SEED, wobble_mutations=True)
        self.complex_name = 'activated_complex'

    @property
    def prescribed_stem_structure(self):
        return (f"({self.b_domain_size_pre_bulge}"
                f".{len(START_CODON)}({self.b_domain_size_post_bulge}"
                f".{self.loop_size}){self.b_domain_size_post_bulge}.{len(START_CODON)}"
                f"){self.b_domain_size_pre_bulge}")

    def get_switch(self, trigger_sequence: str, healthy_sequence: str = None,
                   is_watson_crick: bool = True) -> Tuple[str, float, float, str]:

        if len(trigger_sequence) != self.trigger_length:
            raise ValueError(f"The trigger sequence must be exactly {self.trigger_length} nucleotides long for the prokaryotic switch, but you provided {len(trigger_sequence)} nucleotides.")

        if not set(trigger_sequence).issubset(self.VALID_NUCLEOTIDES):
            invalid_chars = set(trigger_sequence) - self.VALID_NUCLEOTIDES
            raise ValueError(f"The trigger sequence contains invalid characters: {invalid_chars}. Only 'A', 'C', 'G', and 'U' are allowed.")

        assert len(trigger_sequence) == self.trigger_length

        trigger_domain = Domain(trigger_sequence, name="trigger domain", material='rna')
        trigger = TargetStrand([trigger_domain], name="trigger")
        trigger_structure = "." * len(trigger_sequence)
        trigger_complex = TargetComplex([trigger], trigger_structure, name='trigger complex')

        if is_watson_crick:
            trigger_reverse_complement = str(Seq(trigger_sequence).reverse_complement_rna())
            a_domain = Domain(f"{trigger_reverse_complement[:self.a_domain_size]}", name="a", material="rna")
            b_domain_pre_bulge = Domain(
                f"{trigger_reverse_complement[self.a_domain_size:self.a_domain_size + self.b_domain_size_pre_bulge]}",
                name="b_domain_pre_bulge", material="rna")
            b_domain_bulge = Domain(
                f"{trigger_reverse_complement[self.a_domain_size + self.b_domain_size_pre_bulge:self.a_domain_size + self.b_domain_size_pre_bulge + len(START_CODON)]}",
                name="b_domain_bulge", material="rna")
            b_domain_post_bulge = Domain(
                f"{trigger_reverse_complement[-self.b_domain_size_post_bulge:]}",
                name="b_domain_post_bulge", material="rna")
        else:
            a_domain = Domain(f"N{self.a_domain_size}", name="a", material="rna")
            b_domain_pre_bulge = Domain(f"N{self.b_domain_size_pre_bulge}",
                                        name="b_domain_pre_bulge", material="rna")
            b_domain_bulge = Domain(f"N{len(START_CODON)}", name="b_domain_bulge", material="rna")
            b_domain_post_bulge = Domain(f"N{self.b_domain_size_post_bulge}",
                                         name="b_domain_post_bulge", material="rna")

        switch = TargetStrand([
            self.leader, a_domain, b_domain_pre_bulge, b_domain_bulge, b_domain_post_bulge,
            self.loop, ~b_domain_post_bulge, self.bulge, ~b_domain_pre_bulge, self.linker, self.translated_gene
        ], name="switch")
        stop_codon_pattern = Pattern(STOP_CODONS, scope=[~b_domain_pre_bulge])

        switch_complex = TargetComplex([switch], self.switch_structure, name='switch complex')

        if healthy_sequence:
            healthy_sequence_domain = Domain(healthy_sequence, name="healthy trigger domain", material='rna')
            healthy_trigger = TargetStrand([healthy_sequence_domain], name="healthy trigger")
            healthy_trigger_structure = "." * len(healthy_sequence)
            healthy_trigger_complex = TargetComplex([healthy_trigger], healthy_trigger_structure,
                                                    name='healthy trigger complex')

        activated_structure = self._get_activated_structure(switch)
        activated_complex = TargetComplex([trigger, switch], activated_structure, name=self.complex_name)

        reactants = TargetTube(on_targets={switch_complex: 1e-08, trigger_complex: 1e-08},
                               off_targets=SetSpec(max_size=2, exclude=tuple([activated_complex])), name='reactants')

        if healthy_sequence:
            off_target_reactants = TargetTube(on_targets={switch_complex: 1e-08, healthy_trigger_complex: 1e-08},
                                              off_targets=SetSpec(max_size=2), name='off_target_reactants')

        products = TargetTube(on_targets={activated_complex: 1e-08}, off_targets=SetSpec(max_size=2), name='products')

        if healthy_sequence:
            tubes = [reactants, products, off_target_reactants]
        else:
            tubes = [reactants, products]

        my_design = tube_design(tubes=tubes,
                                hard_constraints=[],
                                soft_constraints=[],
                                model=self.model, options=self.options)

        design_results = my_design.run(trials=1)[0]

        switch_designed_strand = [v for k, v in design_results.to_analysis.strands.items() if
                                  v.name == "switch"].pop()
        ensemble_defect = self._get_ensemble_defect(design_results)
        complex_concentration = self._get_trigger_switch_complex_concentration(design_results, self.complex_name)
        switch_mfe_structure = self._get_switch_mfe_structure(design_results, switch_complex)
        return str(switch_designed_strand), ensemble_defect, complex_concentration, switch_mfe_structure

    def _get_activated_structure(self, switch) -> str:
        return ('(' * self.trigger_length + '+' + '.' * self.leader.nt() + ')' * self.trigger_length +
                '.' * (switch.nt() - self.leader.nt() - self.trigger_length))

    def _get_ensemble_defect(self, design_result) -> float:
        return design_result.defects.ensemble_defect

    def _get_trigger_switch_complex_concentration(self, design_result, complex_name) -> float:
        complex_concentrations = design_result.concentrations.table[
            design_result.concentrations.table.complex_name == complex_name].iloc[0]
        target_concentration = complex_concentrations["target_concentration"]
        predicted_concentration = complex_concentrations["concentration"]
        return predicted_concentration / target_concentration

    def _get_switch_mfe_structure(self, design_result, switch_target_strand) -> str:
        return str(mfe(strands=design_result.to_analysis(switch_target_strand), model=self.model)[0].structure)


# Utility function
def check_forbidden_patterns(sequence):
    """Checks if a sequence contains any forbidden patterns."""
    for pattern in FORBIDDEN_PATTERNS:
        if pattern in sequence:
            return True
    return False


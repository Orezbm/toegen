from nupack import *
import pandas as pd
from joblib import Parallel, delayed
from fuzzysearch import find_near_matches
from typing import List, Tuple
from ViennaRNA import RNA
import time
from .switch_generator import SwitchGenerator, ProkaryoticSwitchGenerator

config.threads = 8


class HomologySwitchGenerator:
    def __init__(self, cell_type: str, reporter_gene: str, transcripts_dict: dict):

        self.cell_type = cell_type
        self.reporter_gene = reporter_gene
        self.transcripts_dict = transcripts_dict

    def generate_switch_with_homology(self, triggers: List[str], transcripts_dict: dict,
                                      n_jobs: int = 4) -> pd.DataFrame:
        # Homology search with parallelization
        s_h = time.time()
        homo_res = Parallel(n_jobs=n_jobs)(
            delayed(self.find_homology)(trigger, transcripts_dict) for trigger in triggers)
        e_h = time.time()
        print(f'Homology search time= {e_h - s_h} seconds, n_triggers = {len(triggers)}, cell= {self.cell_type}')

        # Extract top homology sequences with parallelization
        rrf_ranks = self.extract_top_homology_sequences(homo_res)
        homology_sequences = [ranked_df['sequence'].get(0) for ranked_df in rrf_ranks]

        # Generate switches with parallelization
        s_switch = time.time()
        switch_res = Parallel(n_jobs=n_jobs)(delayed(self.generate_switch)(trigger, top_homology_sequence)
                                             for trigger, top_homology_sequence in zip(triggers, homology_sequences)
                                             )
        e_switch = time.time()
        print(f'Switch generation time= {e_switch - s_switch} seconds')

        # Combine results
        results = pd.DataFrame(switch_res, columns=['switch', 'complex_concentration'])
        results['trigger_window'] = triggers
        return results

    def find_homology(self, trigger: str, genome_data: dict) -> List[dict]:
        """
        Searches for homologous sequences in the genome data.
        """
        genes_sub_sequences = []
        for gene_data_dict in genome_data:
            gene_sub_seqs = {}
            mRNA_seq = gene_data_dict['sequence']
            gene_name = gene_data_dict['gene']
            protein = gene_data_dict['protein']
            seq_match_mapping = self.build_homology_map(trigger, mRNA_seq, gene_name, protein)
            if seq_match_mapping:
                gene_sub_seqs[gene_name] = seq_match_mapping
                genes_sub_sequences.append(gene_sub_seqs)
        return genes_sub_sequences

    def build_homology_map(self, trigger: str, seq: str, gene_name: str, protein_name: str) -> List[dict]:
        """
        Builds a homology map of the trigger against a gene sequence.
        """
        sequence_match_mapping = []
        matches = find_near_matches(trigger, seq, max_insertions=0, max_deletions=0, max_l_dist=4)
        for match_obj in matches:
            locus_start = match_obj.start
            locus_end = match_obj.end
            sub_seq = seq[locus_start:locus_end]
            seq_location = (locus_start, locus_end)
            distance = match_obj.dist
            sequence_match_mapping.append({
                'distance': distance, 'idx': seq_location,
                'sequence': sub_seq, 'gene': gene_name, 'protein': protein_name
            })
        return sequence_match_mapping

    def extract_top_homology_sequences(self, triggers_homology_mapping: List[List[dict]]) -> List[pd.DataFrame]:
        """
        Extracts the top-ranked homologous sequences for each trigger.
        """
        homo_dfs = []
        for trigger_homology in triggers_homology_mapping:
            if trigger_homology:
                trigger_all_matches = []
                for gene_homology in trigger_homology:
                    trigger_all_matches.extend(gene_homology.values())
                matches_df = pd.DataFrame(trigger_all_matches[0])

                mfe_dict = {'homologous_trigger_mfe': []}
                for homos in trigger_all_matches[0]:
                    homo_trigger = homos['sequence']
                    structure, mfe = RNA.fold(homo_trigger)
                    mfe_dict['homologous_trigger_mfe'].append(mfe)
                mfe_df = pd.DataFrame(mfe_dict)

                trig_res_df = pd.concat([matches_df, mfe_df], axis=1)
            else:
                trig_res_df = pd.DataFrame(
                    columns=['distance', 'idx', 'sequence', 'gene', 'protein', 'homologous_trigger_mfe'])
            homo_dfs.append(trig_res_df)

        higher = ['homologous_trigger_mfe']
        lower = ['distance']
        rrf_rank = [self.RRF(trigger_homology_df, higher, lower, index='sequence') for trigger_homology_df in homo_dfs]
        return rrf_rank

    def RRF(self, ranking_df: pd.DataFrame, higher_is_better_cols: List[str], lower_is_better_cols: List[str],
            index: str, k: int = 60) -> pd.DataFrame:
        """
        Rank fusion (RRF) algorithm to rank sequences based on multiple metrics.
        """
        ranking_df = ranking_df.copy().reset_index(drop=True)

        for col in higher_is_better_cols:
            ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=False)
        for col in lower_is_better_cols:
            ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=True)

        ranked_columns = higher_is_better_cols + lower_is_better_cols
        ranked_columns = [f'{col}_rank' for col in ranked_columns]
        ranking_df[f'{index}_RRF'] = ranking_df[ranked_columns].apply(lambda row: sum(1 / (k + rank) for rank in row),
                                                                      axis=1)

        return ranking_df.sort_values(by=f'{index}_RRF', ascending=True)

    def generate_switch(self, trigger: str, homologous_sequence: str) -> Tuple[str, float]:
        """
        Generates a switch for the given trigger and homologous sequence.
        """
        if len(trigger) == 30:  # Prokaryotic trigger
            switch_generator = ProkaryoticSwitchGenerator(translated_gene_sequence=self.reporter_gene)
        elif len(trigger) == 23:  # Eukaryotic trigger
            switch_generator = SwitchGenerator(translated_gene_sequence=self.reporter_gene)
        else:
            raise ValueError(f"Invalid trigger length {len(trigger)}. It should be 23 for eukaryotic or 30 for prokaryotic.")

        switch_designed_strand, ensemble_defect, complex_concentration, switch_mfe_structure = switch_generator.get_switch(
            trigger_sequence=trigger, healthy_sequence=homologous_sequence)
        return switch_designed_strand, complex_concentration

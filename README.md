# ToeGen: Toehold Switch Generator with Homology Ranking

`ToeGen` is a Python package designed to generate toehold switches for gene expression control, including the ability to rank homologous sequences. This is useful for applications like synthetic biology and gene expression analysis.

## Features
- Generate eukaryotic and prokaryotic toehold switches.
- Perform homology searches against custom genome sequences.
- Rank homologous sequences based on mfe (minimum free energy) and distance.

## Installation

### Step 1: Install Required Dependencies

Ensure you have the following dependencies installed in your Python environment:

- **Biopython**
- **pandas**
- **joblib**
- **fuzzysearch**
- **ViennaRNA**

These can be installed using `pip`:

```bash
pip install biopython numpy pandas joblib fuzzysearch viennarna
```

### Step 2: Install Nupack

`nupack` is a special requirement that **cannot** be installed via pip. Please follow the official installation guide [here](https://docs.nupack.org/start).

#### Instructions to Install NUPACK:
1. Visit the official [NUPACK download page](https://nupack.org/downloads).
2. Download the appropriate version for your operating system.
3. Follow the installation instructions on the [Getting Started Guide](https://docs.nupack.org/start).
4. After installation, ensure that `nupack` is properly added to your environment.

### Step 3: Install ToeGen

Once the dependencies are installed, you can install `ToeGen`:

```bash
pip install toegen
```

### Usage

Once installed, you can use ToeGen to generate eukaryotic and prokaryotic toehold switches with homology ranking.

- 
- Use Case: Generating a Eukaryotic Toehold Switch
```bash
from toegen import SwitchGenerator

# Example trigger sequence for eukaryotic systems (23 nucleotides)
trigger_sequence = "AUGGCCUAGCGCUAUGCCCUAGG"

# Example reporter gene sequence
reporter_gene = "AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG"

# Initialize the SwitchGenerator for eukaryotic systems
switch_generator = SwitchGenerator(translated_gene_sequence=reporter_gene)

# Generate the toehold switch
switch_strand, ensemble_defect, concentration, mfe_structure = switch_generator.get_switch(trigger_sequence)

# Display the results
print(f"Eukaryotic Switch Strand: {switch_strand}")
print(f"Ensemble Defect: {ensemble_defect}")
print(f"Complex Concentration: {concentration}")
print(f"MFE Structure: {mfe_structure}")
```

- Use Case: Generating a Eukaryotic Toehold Switch
```bash
from toegen import ProkaryoticSwitchGenerator

# Example trigger sequence for prokaryotic systems (30 nucleotides)
trigger_sequence = "AUGGCUCGAUUGCCCUUCUAGGAUCGUAGC"

# Example reporter gene sequence
reporter_gene = "AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG"

# Initialize the ProkaryoticSwitchGenerator for prokaryotic systems
prokaryotic_switch_generator = ProkaryoticSwitchGenerator(translated_gene_sequence=reporter_gene)

# Generate the toehold switch
switch_strand, ensemble_defect, concentration, mfe_structure = prokaryotic_switch_generator.get_switch(trigger_sequence)

# Display the results
print(f"Prokaryotic Switch Strand: {switch_strand}")
print(f"Ensemble Defect: {ensemble_defect}")
print(f"Complex Concentration: {concentration}")
print(f"MFE Structure: {mfe_structure}")
```

- Use Case: Generating Toehold Switches with Homology Ranking
```bash
from toegen import HomologySwitchGenerator

# Example genome data simulating transcript sequences (replace with real data if available)
transcripts_dict = [
    {
        'gene': 'GeneA',
        'protein': 'ProteinA',
        'sequence': 'AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG'
    },
    {
        'gene': 'GeneB',
        'protein': 'ProteinB',
        'sequence': 'AUGGCUCGAUUGCCCUUCUAGGAUCGUAGCUAGGAUC'
    }
]

# Example list of trigger sequences
triggers = [
    "AUGGCCUAGCGCUAUGCCCUAGG",  # Eukaryotic trigger (23 nucleotides)
    "AUGGCUCGAUUGCCCUUCUAGGAUCGUAGC"  # Prokaryotic trigger (30 nucleotides)
]

# Example reporter gene sequence
reporter_gene = "AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG"

# Initialize the HomologySwitchGenerator class for Homo sapiens (eukaryotic)
homology_switch_generator = HomologySwitchGenerator(
    cell_type="Homo sapiens",  # Specify the cell type
    reporter_gene=reporter_gene,  # Specify the reporter gene sequence
    transcripts_dict=transcripts_dict  # Pass in mock genome data
)

# Generate toehold switches with homology search and ranking
results = homology_switch_generator.generate_switch_with_homology(triggers, homology_switch_generator.transcripts_dict)

# Display the results (ranked homology matches and switch generation)
print(results)
```

### License

This project is licensed under the terms of the MIT license.
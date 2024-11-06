from genedesign.seq_utils.Translate import Translate
from genedesign.seq_utils.hairpin_counter import hairpin_counter
from genedesign.seq_utils.calc_edit_distance import calculate_edit_distance
from genedesign.models.rbs_option import RBSOption
from typing import Set

import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import random


class RBSChooser:
    """
    A class to choose the best RBS for a given CDS sequence based on multiple criteria:
    1. Excluding RBS options from an 'ignores' list.
    2. Reducing secondary structure occlusion between the RBS and the input CDS.
    3. Maximizing the similarity of the input peptide to the RBS source gene's peptide.

    The class uses a multi-objective optimization strategy to select the most suitable RBS option.

    Attributes:
        rbs_options (list[RBSOption]): A list of available RBSOption instances.
        translator (Translate): A Translate instance to convert DNA sequences to protein sequences.
    """

    def __init__(self):
        """
        Initializes the RBSChooser instance and prepares the list of RBSOptions.
        """
        self.rbs_options = []  # List of RBSOptions
        self.translator = Translate()  # Translator to handle DNA to protein conversion
        self.translator.initiate()
        genbank_file = "genedesign/data/sequence.gb"
        proteomics_file = "genedesign/data/511145-WHOLE_ORGANISM-integrated.txt"
        top_percentage = 5
        # Automatically generate rbs_cds_dict as part of initialization
        genes_info = self.extract_genes_info(genbank_file)
        pruned_proteomics_data = self.prune_proteomics_data(proteomics_file, top_percentage)
        self.rbs_cds_dict = self.create_rbs_cds_dict(pruned_proteomics_data, genes_info)

        # Use the refined rbs_cds_dict to populate RBSOptions
        self.initiate()

    def extract_genes_info(self, genbank_file):
        """
        Extracts gene, UTR, and CDS information from the GenBank file.
        """
        gene_dict = defaultdict(dict)
        for record in SeqIO.parse(genbank_file, "genbank"):
            for feature in record.features:
                if feature.type == "gene":
                    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                    gene_name = feature.qualifiers.get("gene", [None])[0]

                    cds_feature = None
                    for cds in record.features:
                        if cds.type == "CDS" and cds.qualifiers.get("locus_tag") == [locus_tag]:
                            cds_feature = cds
                            break

                    if cds_feature:
                        start, end = cds_feature.location.start, cds_feature.location.end
                        strand = cds_feature.location.strand
                        if strand == 1:
                            utr_start = max(0, start - 50)
                            utr_seq = record.seq[utr_start:start]
                        else:
                            utr_start = end
                            utr_seq = record.seq[utr_start:utr_start + 50].reverse_complement()

                        cds_seq = cds_feature.extract(record.seq)
                        gene_dict[locus_tag] = {
                            "gene": gene_name,
                            "UTR": utr_seq,
                            "CDS": cds_seq
                        }
        return gene_dict

    def prune_proteomics_data(self, file_path, top_percentage=5):
        """
        Prunes the proteomics data to the top specified percentage.
        """
        proteomics_data = pd.read_csv(file_path, sep="\t", skiprows=11)
        if 'abundance' in proteomics_data.columns:
            sorted_data = proteomics_data.sort_values(by='abundance', ascending=False)
        else:
            raise KeyError("The 'abundance' column was not found in the dataset.")

        top_rows = int(len(sorted_data) * (top_percentage / 100.0))
        pruned_data = sorted_data.head(top_rows)
        return pruned_data

    def create_rbs_cds_dict(self, pruned_data, genes_info):
        """
        Creates a dictionary by merging pruned proteomics data with genes_info.
        """
        merged_dict = {}
        for _, row in pruned_data.iterrows():
            locus_tag = row['#string_external_id'].replace('511145.', '')
            if locus_tag in genes_info:
                merged_dict[locus_tag] = genes_info[locus_tag]
        return merged_dict

    def initiate(self):
        """
        Initializes the RBSChooser with RBSOptions based on the filtered rbs_cds_dict.
        """
        def calculate_first_six_aas(cds: str) -> str:
            if len(cds) < 18:
                raise ValueError("CDS length must be at least 18 nucleotides to calculate first six amino acids.")
            return self.translator.run(cds[:18])

        for locus_tag, gene_data in self.rbs_cds_dict.items():
            utr = str(gene_data['UTR'])
            cds = str(gene_data['CDS'])
            gene_name = gene_data['gene']

            if '...' in utr or '...' in cds:
                continue

            first_six_aas = calculate_first_six_aas(cds)
            rbs_option = RBSOption(
                utr=utr,
                cds=cds,
                gene_name=gene_name,
                first_six_aas=first_six_aas
            )
            self.rbs_options.append(rbs_option)

    def run(self, cds: str, ignores: Set[RBSOption], num_samples=100, max_iterations=100) -> RBSOption:
        """
        Uses a Monte Carlo sampling approach to select the best RBSOption for a given CDS sequence.
        
        Parameters:
            cds (str): The coding sequence for which an RBS needs to be selected.
            ignores (Set[RBSOption]): A set of RBSOptions to be ignored in the selection process.
            num_samples (int): Number of random samples to evaluate per iteration.
            max_iterations (int): Maximum number of Monte Carlo sampling iterations.
        
        Returns:
            RBSOption: The selected RBSOption that best fits the criteria, or None if no valid options are available.
        """
        # Validate CDS input
        if not isinstance(cds, str) or len(cds) < 18 or len(cds) % 3 != 0:
            raise ValueError("CDS sequence must be a string of length at least 18 and a multiple of 3.")

        # Step 1: Filter out ignored RBS options
        available_rbs_options = [rbs for rbs in self.rbs_options if rbs not in ignores]
        if not available_rbs_options:
            return None

        # Track the best RBS option and score found
        best_rbs_option = None
        best_score = float('inf')
        input_peptide = self.translator.run(cds[:18])  # First 6 codons of the input CDS

        for _ in range(max_iterations):
            # Step 2: Randomly sample a subset of RBS options
            sampled_rbs_options = random.sample(available_rbs_options, min(num_samples, len(available_rbs_options)))
            
            # Evaluate each sampled RBS option
            for rbs in sampled_rbs_options:
                # Calculate hairpin count and edit distance without caching
                hairpin_count, _ = hairpin_counter(rbs.utr + cds)
                edit_distance = calculate_edit_distance(input_peptide, rbs.first_six_aas)

                # Calculate score (lower is better)
                total_score = hairpin_count + edit_distance

                # Update best option if this one is better
                if total_score < best_score:
                    best_rbs_option = rbs
                    best_score = total_score

            # Early stopping condition: if best score stabilizes below a threshold
            if best_score < 5:  # This threshold is arbitrary and can be tuned
                break

        return best_rbs_option

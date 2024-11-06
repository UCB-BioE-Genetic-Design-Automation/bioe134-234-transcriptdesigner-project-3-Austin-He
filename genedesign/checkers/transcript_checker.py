from genedesign.seq_utils.reverse_complement import reverse_complement

class TranscriptChecker:
    """
    A class to check for the presence of specific transcript motifs in a DNA sequence.

    This class uses domain knowledge to identify common features associated with transcript motifs,
    such as the presence of initiation codons, specific binding sites, and other sequences known
    to be involved in transcript processing and recognition.

    Attributes:
        initiation_codons: A list of initiation codons typically found in transcripts.
        polyadenylation_signal: A common polyadenylation signal sequence found in transcripts.
        splice_sites: A list of common splice donor and acceptor sites.
        promoter_regions: A list of common promoter motifs associated with transcript initiation.
        terminator_sequences: A list of sequences often found at the end of transcription units.
        ribosome_binding_sites: A list of common ribosome binding site sequences (Shine-Dalgarno sequences).
        enhancers: A list of common enhancer sequences that aid in transcription initiation.
        intron_motifs: A list of sequences typically found at the start and end of introns.
    """

    def __init__(self):
        """
        Initializes the TranscriptChecker with known domain-specific sequences.
        """
        self.initiation_codons = ["ATG"]  # Common initiation codon for transcripts
        self.polyadenylation_signal = "AATAAA"  # Common polyadenylation signal sequence
        self.splice_sites = {
            "donor": "GT",  # Common splice donor site
            "acceptor": "AG"  # Common splice acceptor site
        }
        self.promoter_regions = ["TATAAT", "TTGACA"]  # Common promoter motifs (TATA box, -35 region)
        self.terminator_sequences = ["TTATT", "TATGTT"]  # Common terminator sequences
        self.ribosome_binding_sites = ["AGGAGG"]  # Shine-Dalgarno sequence for ribosome binding
        self.enhancers = ["CAAT", "GGGCGG"]  # Common enhancer sequences
        self.intron_motifs = {
            "start": "GT",  # Intron start sequence
            "end": "AG"  # Intron end sequence
        }

    def run(self, seq):
        """
        Checks if the given DNA sequence contains transcript-related features.

        The method looks for the presence of initiation codons, polyadenylation signals, splice sites,
        promoter regions, terminator sequences, ribosome binding sites, enhancer sequences, and intron motifs.
        It also evaluates both the input sequence and its reverse complement.

        Parameters:
            seq (str): A DNA sequence to check.

        Returns:
            tuple: (bool, str or None)
                - bool: True if no transcript motif is found, False if a motif is found.
                - str: The motif sequence if found, None otherwise.
        """
        # Convert the sequence to uppercase and compute its reverse complement.
        seq = seq.upper()
        rc = reverse_complement(seq)
        combined = [seq, rc]  # List containing the original sequence and its reverse complement.

        # Check for transcript features in both original and reverse complement sequences.
        for strand in combined:
            # Check for initiation codons
            for codon in self.initiation_codons:
                if codon in strand:
                    start_index = strand.find(codon)
                    return False, strand[start_index:start_index + 15]  # Return the sequence around the initiation codon

            # Check for polyadenylation signal
            if self.polyadenylation_signal in strand:
                signal_index = strand.find(self.polyadenylation_signal)
                return False, strand[signal_index:signal_index + len(self.polyadenylation_signal)]  # Return the polyadenylation signal

            # Check for splice sites (donor and acceptor)
            if self.splice_sites["donor"] in strand and self.splice_sites["acceptor"] in strand:
                donor_index = strand.find(self.splice_sites["donor"])
                acceptor_index = strand.find(self.splice_sites["acceptor"], donor_index)
                if acceptor_index != -1:
                    return False, strand[donor_index:acceptor_index + len(self.splice_sites["acceptor"])]  # Return splice site sequence

            # Check for promoter regions
            for promoter in self.promoter_regions:
                if promoter in strand:
                    promoter_index = strand.find(promoter)
                    return False, strand[promoter_index:promoter_index + len(promoter)]  # Return the promoter sequence

            # Check for terminator sequences
            for terminator in self.terminator_sequences:
                if terminator in strand:
                    terminator_index = strand.find(terminator)
                    return False, strand[terminator_index:terminator_index + len(terminator)]  # Return the terminator sequence

            # Check for ribosome binding sites (Shine-Dalgarno sequences)
            for rbs in self.ribosome_binding_sites:
                if rbs in strand:
                    rbs_index = strand.find(rbs)
                    return False, strand[rbs_index:rbs_index + len(rbs)]  # Return the ribosome binding site sequence

            # Check for enhancer sequences
            for enhancer in self.enhancers:
                if enhancer in strand:
                    enhancer_index = strand.find(enhancer)
                    return False, strand[enhancer_index:enhancer_index + len(enhancer)]  # Return the enhancer sequence

            # Check for intron motifs (start and end)
            if self.intron_motifs["start"] in strand and self.intron_motifs["end"] in strand:
                intron_start_index = strand.find(self.intron_motifs["start"])
                intron_end_index = strand.find(self.intron_motifs["end"], intron_start_index)
                if intron_end_index != -1:
                    return False, strand[intron_start_index:intron_end_index + len(self.intron_motifs["end"])]  # Return intron sequence

        return True, None  # No transcript motif detected in the sequence


if __name__ == "__main__":
    checker = TranscriptChecker()

    transcript = "ATGCGTACGTAGCTG"
    transcriptBroken = "GTACCGTATCAGGTC"

    # Example checks
    result, motif = checker.run(transcript)
    print(f"transcript: {result}, Motif: {motif}")  # Expected output: False (contains motif)

    result, motif = checker.run(transcriptBroken)
    print(f"transcriptBroken: {result}, Motif: {motif}")  # Expected output: True (no motif)

    # Example of testing a list of sequences
    seq_list = [
        "ATGCGTACGTAGCTG",
        "GTACCGTATCAGGTC",
        "ATGCGTACGTACGTA",
        "TACGTACGTAGCGTT",
        "GGGCGTACGTAGGGG"
    ]

    print("\n>> Testing transcript motifs (mixed results)")
    for seq in seq_list:
        result, motif = checker.run(seq)
        print(f"Result: {result}, Motif: {motif}")  # Outputs either True (no motif) or False (contains motif)

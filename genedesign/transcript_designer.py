from genedesign.rbs_chooser import RBSChooser 
from genedesign.models.transcript import Transcript
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
import random

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using a Monte Carlo approach for each amino acid.
    """

    def __init__(self):
        self.aminoAcidToCodons = {}
        self.rbsChooser = None
        self.codon_checker = None
        self.forbidden_sequence_checker = None
        self.internal_promoter_checker = None

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # Initialize the checkers
        self.codon_checker = CodonChecker()
        self.codon_checker.initiate()
        self.forbidden_sequence_checker = ForbiddenSequenceChecker()
        self.internal_promoter_checker = PromoterChecker()
        self.internal_promoter_checker.initiate()  # Ensure promoter checker is properly initiated

        # Codon table with possible codons for each amino acid (for E. coli)
        self.aminoAcidToCodons = {
            'A': ["GCT", "GCC", "GCA", "GCG"], 'C': ["TGT", "TGC"], 'D': ["GAT", "GAC"],
            'E': ["GAA", "GAG"], 'F': ["TTT", "TTC"], 'G': ["GGT", "GGC", "GGA", "GGG"],
            'H': ["CAT", "CAC"], 'I': ["ATT", "ATC", "ATA"], 'K': ["AAA", "AAG"],
            'L': ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], 'M': ["ATG"],
            'N': ["AAT", "AAC"], 'P': ["CCT", "CCC", "CCA", "CCG"],
            'Q': ["CAA", "CAG"], 'R': ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
            'S': ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], 'T': ["ACT", "ACC", "ACA", "ACG"],
            'V': ["GTT", "GTC", "GTA", "GTG"], 'W': ["TGG"], 'Y': ["TAT", "TAC"]
        }

    def run(self, peptide: str, ignores: set = None) -> Transcript:
        """
        Optimizes the transcript using a Monte Carlo approach and then applies a sliding window adjustment.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The optimized transcript object, if found.
        """
        if ignores is None:
            ignores = set()

        # Run Monte Carlo to find the best possible transcript design
        optimized_transcript = self.monte_carlo_transcript_optimization(peptide, ignores=ignores)

        # If a transcript is found, return it; otherwise, print a failure message
        if optimized_transcript:
            print(f"Optimized transcript found: {optimized_transcript}")
            return optimized_transcript
        else:
            print("No valid transcript found after Monte Carlo optimization.")
            return None

    def check_transcript_with_fix(self, transcript) -> bool:
        if not transcript:
            return False

        dna_sequence = ''.join(transcript.codons)
        all_checks_passed = True

        # Forbidden Sequence Check
        forbidden_valid, forbidden_seq = self.forbidden_sequence_checker.run(dna_sequence)
        if not forbidden_valid:
            print(f"Forbidden sequence found: {forbidden_seq}")
            transcript = self.fix_forbidden_sequence(transcript, forbidden_seq)
            if not transcript:
                return False
            dna_sequence = ''.join(transcript.codons)
            all_checks_passed = False
        print(f"Forbidden sequence check: {forbidden_valid}")

        # Hairpin Check
        hairpin_valid, hairpin_structure = hairpin_checker(dna_sequence)
        if not hairpin_valid:
            print(f"Hairpin structure found: {hairpin_structure}")
            transcript = self.fix_hairpin(transcript, hairpin_structure)
            if not transcript:
                return False
            dna_sequence = ''.join(transcript.codons)
            all_checks_passed = False
        print(f"Hairpin check: {hairpin_valid}")

        # Promoter Check
        if self.internal_promoter_checker is None or self.internal_promoter_checker.pwm is None:
            print("Error: PromoterChecker is not properly initialized or PWM is missing.")
            return False
        print(f"Internal promoter check: {self.internal_promoter_checker}")
        promoter_valid, promoter_seq = self.internal_promoter_checker.run(dna_sequence)
        while not promoter_valid:
            print(f"Internal promoter sequence found: {promoter_seq}")
            transcript = self.fix_promoter_sequence(transcript, promoter_seq)
            if not transcript:
                return False
            dna_sequence = ''.join(transcript.codons)
            promoter_valid, promoter_seq = self.internal_promoter_checker.run(dna_sequence)
        print(f"Promoter check after fixing: {promoter_valid}")

        # Codon Check
        codon_valid, diversity, rare_count, cai = self.codon_checker.run(transcript.codons)
        print(f"Codon check: Diversity={diversity}, Rare Codons={rare_count}, CAI={cai}")
        if not codon_valid:
            print(f"Codon check failed: Diversity={diversity}, Rare Codons={rare_count}, CAI={cai}")
            all_checks_passed = False

        # Return whether all checks passed
        return all_checks_passed

    def monte_carlo_transcript_optimization(self, peptide: str, ignores: set, max_iterations=1000):
        """
        Runs a Monte Carlo simulation to optimize the transcript design.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
            max_iterations (int): The maximum number of iterations to attempt.
        
        Returns:
            Transcript: The optimized transcript object, if found.
        """
        best_transcript = None
        best_score = float('-inf')

        for iteration in range(max_iterations):
            print(f"Iteration {iteration + 1} of {max_iterations}")
            transcript = self.generate_transcript(peptide, ignores)

            if transcript:
                # Evaluate the transcript
                _, diversity, rare_count, cai = self.codon_checker.run(transcript.codons)
                score = (0.5 * diversity) + (0.3 * cai) - (0.2 * rare_count)
                print(f"Transcript score: {score}")

                # If transcript doesn't meet thresholds, try to adjust
                if diversity < 0.5 or rare_count > 3 or cai < 0.2:
                    transcript = self.sliding_window_adjustment(transcript, window_size=3)
                    if transcript:
                        _, diversity, rare_count, cai = self.codon_checker.run(transcript.codons)
                        score = (0.5 * diversity) + (0.3 * cai) - (0.2 * rare_count)
                        print(f"Adjusted Transcript score: {score}")

                # Keep track of the best transcript
                if score > best_score:
                    best_transcript = transcript
                    best_score = score

                # If we have a perfect score, stop early
                if best_score >= 1.0:
                    break

        if best_transcript:
            print(f"Best transcript found with score {best_score}")
        else:
            print("Max iterations reached without finding a valid transcript.")
        return best_transcript

    def generate_transcript(self, peptide: str, ignores: set) -> Transcript:
        """
        Generates a transcript from the given peptide sequence.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.

        Returns:
            Transcript: A generated transcript object.
        """
        # Translate peptide to codons
        codons = [random.choice(self.aminoAcidToCodons[aa]) for aa in peptide]

        # Append the stop codon (TAA in this case)
        codons.append("TAA")

        # Build the CDS from the codons
        cds = ''.join(codons)

        # Choose an RBS, pass ignores set to ensure proper filtering
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Create a Transcript object
        transcript = Transcript(selectedRBS, peptide, codons)

        # Run the checking logic
        if self.check_transcript_with_fix(transcript):
            return transcript
        else:
            return None

    def sliding_window_adjustment(self, transcript: Transcript, window_size=9):
        """
        Adjusts the transcript using a sliding window approach to optimize it.
        
        Parameters:
            transcript (Transcript): The initial transcript to adjust.
            window_size (int): The size of the sliding window.
        
        Returns:
            Transcript: The adjusted transcript object.
        """
        if not transcript:
            return None

        codons = transcript.codons[:]
        for i in range(0, len(codons) - window_size + 1):
            window = codons[i:i + window_size]
            for j in range(len(window)):
                aa = self.translator.run(window[j])
                if aa in self.aminoAcidToCodons:
                    new_codon = random.choice(self.aminoAcidToCodons[aa])
                    print(f"Sliding window adjustment: Replacing codon {window[j]} at index {i + j} with {new_codon}")
                    codons[i + j] = new_codon

        # Create a new transcript with adjusted codons
        adjusted_transcript = Transcript(transcript.rbs, transcript.peptide, codons)
        return adjusted_transcript

    def fix_forbidden_sequence(self, transcript: Transcript, forbidden_seq: str):
        """
        Fixes the forbidden sequence in the transcript by replacing the codons involved using a Monte Carlo and sliding window approach.
        
        Parameters:
            transcript (Transcript): The transcript to fix.
            forbidden_seq (str): The forbidden sequence to be replaced.
        """
        if not transcript:
            return None

        dna_sequence = ''.join(transcript.codons)
        start_index = dna_sequence.find(forbidden_seq)
        if start_index != -1:
            end_index = start_index + len(forbidden_seq)
            window_size = 9
            for i in range(start_index, end_index, 3):
                for _ in range(10):  # Monte Carlo iterations to find a suitable replacement
                    window_start = max(0, i - window_size // 2)
                    window_end = min(len(dna_sequence), i + window_size // 2)
                    window = dna_sequence[window_start:window_end]
                    aa = self.translator.run(dna_sequence[i:i+3])
                    if aa in self.aminoAcidToCodons:
                        new_codon = random.choice(self.aminoAcidToCodons[aa])
                        transcript.codons[i // 3] = new_codon
                        print(f"Fixed forbidden sequence: Replacing codon at index {i // 3} with {new_codon}")
                        break
        return transcript

    def fix_hairpin(self, transcript: Transcript, hairpin_structure: str):
        """
        Fixes the hairpin structure in the transcript by replacing the codons involved using a Monte Carlo and sliding window approach.
        
        Parameters:
            transcript (Transcript): The transcript to fix.
            hairpin_structure (str): The hairpin structure sequence to be replaced.
        """
        if not transcript:
            return None

        dna_sequence = ''.join(transcript.codons)
        start_index = dna_sequence.find(hairpin_structure)
        if start_index != -1:
            end_index = start_index + len(hairpin_structure)
            window_size = 9
            for i in range(start_index, end_index, 3):
                for _ in range(10):  # Monte Carlo iterations to find a suitable replacement
                    window_start = max(0, i - window_size // 2)
                    window_end = min(len(dna_sequence), i + window_size // 2)
                    window = dna_sequence[window_start:window_end]
                    aa = self.translator.run(dna_sequence[i:i+3])
                    if aa in self.aminoAcidToCodons:
                        new_codon = random.choice(self.aminoAcidToCodons[aa])
                        transcript.codons[i // 3] = new_codon
                        print(f"Fixed hairpin structure: Replacing codon at index {i // 3} with {new_codon}")
                        break
        return transcript

    def fix_promoter_sequence(self, transcript: Transcript, promoter_seq: str):
        """
        Fixes the internal promoter sequence in the transcript by replacing the codons involved using a Monte Carlo and sliding window approach.
        
        Parameters:
            transcript (Transcript): The transcript to fix.
            promoter_seq (str): The promoter sequence to be replaced.
        """
        if not transcript:
            return None

        dna_sequence = ''.join(transcript.codons)
        start_index = dna_sequence.find(promoter_seq)
        if start_index != -1:
            end_index = start_index + len(promoter_seq)
            window_size = 9
            while not self.internal_promoter_checker.run(dna_sequence)[0]:
                for i in range(start_index, end_index, 3):
                    for _ in range(10):  # Monte Carlo iterations to find a suitable replacement
                        window_start = max(0, i - window_size // 2)
                        window_end = min(len(dna_sequence), i + window_size // 2)
                        window = dna_sequence[window_start:window_end]
                        aa = self.translator.run(dna_sequence[i:i+3])
                        if aa in self.aminoAcidToCodons:
                            new_codon = random.choice(self.aminoAcidToCodons[aa])
                            transcript.codons[i // 3] = new_codon
                            dna_sequence = ''.join(transcript.codons)
                            print(f"Fixed promoter sequence: Replacing codon at index {i // 3} with {new_codon}")
                            break
                promoter_valid, promoter_seq = self.internal_promoter_checker.run(dna_sequence)
                if promoter_valid:
                    break
        return transcript
                
if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    
    peptide = "MYPFIRTARMTV"
    
    #this only has 13 characters and the diversity calculation divides num of unique codons by 62
    # so its impossible to meet the 0.5 threshold with this length peptide because there can't be enough codons
    # but we can get it to have 0 rare codons easily.
    
    #peptide = 'MRRGLVIVGHGSQLNHYREVMELHRKRIEESGAFDEVKIAFAARKRRPMPDEAIREMNCDIIYVVPLFISYGLHVTEDLPDLLGFPRGRGIKEGEFEGKKVVICEPIGEDYFVTYAILNSVFRIGRDGKGEE'
    
    #this has a hundred something characters and gets a higher diversity at aroud 0.7 to 0.8 because we have space to use more codons
    #its also way more likely to have more rare codons because its bigger, it gets around 10
    
    #I feel like the diversity metric has to be normalized somehow and dividng by 62 isn't right. t
    
    #sometimes when trying to fix an issue with hairpin or promoter we get returned None. 
    # this causes AttributeError: 'NoneType' object has no attribute 'codons'
    # we could add some logic around to stop this from causing an attribute error.
     #sometimes it also gets stuck on internal promoters and keeps trying to fix them forever. 
    #although in some earlier version this didn't happen
    
    designer = TranscriptDesigner()
    designer.initiate()

    # Run Monte Carlo with Sliding Window optimization
    optimized_transcript = designer.run(peptide)
    
    # Print out the optimized transcript information
    if optimized_transcript:
        print(optimized_transcript)
    else:
        print("No valid transcript found.")
    
    # Print internal state history
    print("\nInternal State History:")
    for state in designer.history:
        print(state)

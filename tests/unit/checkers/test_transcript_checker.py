import pytest
from genedesign.seq_utils.reverse_complement import reverse_complement
from genedesign.checkers.transcript_checker import TranscriptChecker

@pytest.fixture
def transcript_checker():
    return TranscriptChecker()

def test_initiation_codon(transcript_checker):
    seqs_with_initiation_codon = [
        "ATGCGTACGTAGCTG",
        "GGGATGCGGATAGGAG",
        "TTGATGCTCAGTGTATG",
    ]
    print("\n>> Testing sequences with initiation codon (expected False)")
    for seq in seqs_with_initiation_codon:
        result, motif = transcript_checker.run(seq)
        print(f"Result: {result} on {seq}, Motif: {motif}")
        assert result == False  # Expecting sequences with initiation codon to return False

def test_no_transcript_motif(transcript_checker):
    random_seqs = [
        "AAACTGTAATCCACCACAAGTCAAGCCAT",
        "GCCTCTCTGAGACGCCGTATGAATTAATA",
        "GTAAACTTTGCGCGGGTTCACTGCGATCC",
        "TTCAGTCTCGTCCAAGGGCACAATCGAAT",
        "ATCCCCCGAAGTTTAGCAGGTCGTGAGGT",
        "TCATGGAGGCTCTCGTTCATCCCGTGGGA",
    ]
    print("\n>> Testing sequences with no expected transcript motifs (expected True)")
    for seq in random_seqs:
        result, motif = transcript_checker.run(seq)
        print(f"Result: {result} on {seq}, Motif: {motif}")
        assert result == True  # Expecting no motifs to return True

def test_polyadenylation_signal(transcript_checker):
    seq_with_polyA_signal = [
        "AAATAATAAAGGCTTAG",
        "CCGTAATAAATAGCTGC",
        "GGGTAATAAATGCCATA",
    ]
    print("\n>> Testing sequences with polyadenylation signals (expected False)")
    for seq in seq_with_polyA_signal:
        result, motif = transcript_checker.run(seq)
        print(f"Result: {result} on {seq}, Motif: {motif}")
        assert result == False  # Expecting sequences with polyadenylation signal to return False

def test_ribosome_binding_sites(transcript_checker):
    seqs_with_rbs = [
        "GGAAGGAGGTAGCTGA",
        "TTTAGGAGGATCGACT",
        "AGGAGGCGTGAGTAGT",
    ]
    print("\n>> Testing sequences with ribosome binding sites (expected False)")
    for seq in seqs_with_rbs:
        result, motif = transcript_checker.run(seq)
        print(f"Result: {result} on {seq}, Motif: {motif}")
        assert result == False  # Expecting sequences with RBS to return False

def test_splice_sites(transcript_checker):
    seqs_with_splice_sites = [
        "GTAGCTGCTAGAGCTAG",
        "GTGGGATCAGAGGTAAG",
        "GTACTGAGAGGTCAGAG",
    ]
    print("\n>> Testing sequences with splice sites (expected False)")
    for seq in seqs_with_splice_sites:
        result, motif = transcript_checker.run(seq)
        print(f"Result: {result} on {seq}, Motif: {motif}")
        assert result == False  # Expecting sequences with splice sites to return False

def test_enhancer_sequences(transcript_checker):
    seqs_with_enhancers = [
        "CAATGCGTGCTAGCTA",
        "GGGCGGTAGCTAGGTC",
        "TACCAATGGCCGGGCGG",
    ]
    print("\n>> Testing sequences with enhancer sequences (expected False)")
    for seq in seqs_with_enhancers:
        result, motif = transcript_checker.run(seq)
        print(f"Result: {result} on {seq}, Motif: {motif}")
        assert result == False  # Expecting sequences with enhancer to return False

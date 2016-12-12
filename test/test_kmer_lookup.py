from bioseq.kmer_lookup import KmerLookup
from nose.tools import eq_

def test_kmer_lookup():
    sequences = {"a": "xxxABCyyy", "b": "!!ABC$$"}
    kmer_lookup = KmerLookup(sequence_dictionary=sequences, min_kmer_size=3)
    kmer_lookup.index()
    eq_(kmer_lookup.find_occurrences("ABC"),
        {"a": {3}, "b": {2}})

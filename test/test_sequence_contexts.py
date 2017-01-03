from nose.tools import eq_
from bioseq.sequence_contexts import SequenceContextGenerator, SequenceContext

def test_sequence_context_generator():
    sequences = {"a": "xxxABCyyy", "b": "!!ABC$$"}
    gen = SequenceContextGenerator(
        sequence_dictionary=sequences,
        min_kmer_size=3,
        n_upstream=4,
        n_downstream=3,
        pad_context=True)
    gen.index()
    query_sequences = ["ABC", "DOES_NOT_OCCUR"]
    expected = {
        "ABC": [
            SequenceContext(
                sequence_id="a",
                offset=3,
                upstream="-xxx",
                downstream="yyy"),
            SequenceContext(
                sequence_id="b",
                offset=2,
                upstream="--!!",
                downstream="$$-")
        ],
        "DOES_NOT_OCCUR": []
    }
    eq_(
        gen.sequence_contexts_dictionary(
            query_sequences,
            sort_results=True), expected)

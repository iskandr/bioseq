# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import
from collections import namedtuple

from .kmer_lookup import KmerLookup

SequenceContext = namedtuple(
    "SequenceContext",
    [
        "sequence_id",
        "offset",
        "upstream",
        "downstream"
    ])

class SequenceContextGenerator(KmerLookup):
    """
    In addition looking up where subsequences occur in a group of longer
    sequences, we sometimes also want to know the surrounding sequence
    context. This class lookups up the sequences before/after query hits
    and returns them as SequenceContext namedtuples.
    """
    def __init__(
            self,
            sequence_dictionary,
            min_kmer_size,
            n_upstream,
            n_downstream,
            pad_context=True):
        KmerLookup.__init__(self, sequence_dictionary, min_kmer_size)
        self.n_upstream = n_upstream
        self.n_downstream = n_downstream
        self.pad_context = pad_context

    def sequence_contexts(self, query):
        """
        Generator which yields SequenceContext objects.
        """
        original_sequence_dictionary = self.sequence_dictionary
        n_upstream = self.n_upstream
        n_downstream = self.n_downstream
        pad_context = self.pad_context
        id_to_offsets = self.find_occurrences(query)
        for (identifier, offsets) in id_to_offsets.items():
            for offset in offsets:
                full_sequence = original_sequence_dictionary[identifier]
                upstream = full_sequence[:offset]
                downstream = full_sequence[offset + len(query):]
                if pad_context:
                    upstream = '-' * (n_upstream - len(upstream)) + upstream
                    downstream = downstream + '-' * (n_downstream - len(downstream))
                yield SequenceContext(
                    sequence_id=identifier,
                    offset=offset,
                    upstream=upstream,
                    downstream=downstream)

    def sequence_contexts_dictionary(
            self,
            queries,
            sort_results=False):
        """
        Returns dictionary mapping each query to a list of
        SequenceContext objects.
        """
        query_to_sequence_contexts = {}
        for query in queries:
            results = list(self.sequence_contexts(query))
            if sort_results:
                results.sort(key=lambda x: (x.sequence_id, x.offset))
            query_to_sequence_contexts[query] = results
        return query_to_sequence_contexts

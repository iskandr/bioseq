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
from collections import defaultdict

from six.moves import range

def _list_to_dict(list_of_pairs):
    """
    Given a list of pairs (id, offset) return a dictionary from
    identifiers to set of offsets.
    """
    result_dict = defaultdict(set)
    for identifier, offset in list_of_pairs:
        result_dict[identifier].add(offset)
    return result_dict

class KmerIndex(object):
    def __init__(self, sequence_dictionary, k):
        """
        Given a dictionary from sequence IDs to sequences,
        returns a dictionary from kmers to pairs of (sequence ID, position).
        """
        self.k = k
        self.kmer_index_dict = self._create_index(sequence_dictionary, k)

    def _create_index(self, sequence_dictionary, k):
        kmers = defaultdict(list)
        for (key, sequence) in sequence_dictionary.items():
            for i in range(len(sequence) - k):
                kmer = sequence[i:i + k]
                kmers[kmer].append((key, i))
        return kmers

    def find_subsequence_occurrences(self, query):
        """
        Returns a dictionary of sequence identifiers mapping to sets of
        offsets where the query subsequence occurs.
        """
        k = self.k
        kmer_index_dict = self.kmer_index_dict
        if len(query) < k:
            raise ValueError(
                "Query sequence '%s' must be at least %d characters" % (
                    query, k))

        hits = _list_to_dict(kmer_index_dict.get(query[:k], []))
        for i in range(1, len(query) - k):
            new_hits = defaultdict(set)
            for identifier, offset in kmer_index_dict.get(query[i:i + k], []):
                if identifier in hits and offset - 1 in hits[identifier]:
                    new_hits[identifier].add(offset)
            if len(new_hits) == 0:
                break
            hits = new_hits
        return {
            identifier: {offset - len(query) + k for offset in offsets}
            for (identifier, offsets) in hits.items()
        }

    def find_contexts(
            self,
            queries,
            original_sequence_dictionary,
            n_before,
            n_after,
            padding=True):
        query_to_id_to_contexts = {}
        for query in queries:
            id_to_offsets = self.find_subsequence_occurrences(query)
            id_to_contexts = defaultdict(list)
            for (identifier, offset) in id_to_offsets:
                full_sequence = original_sequence_dictionary[identifier]
                before = full_sequence[:offset]
                after = full_sequence[offset + len(query):]
                if padding:
                    before = '-' * (n_before - len(before)) + before
                    after = after + '-' * (n_after - len(after))
                id_to_contexts[identifier].append((before, after))
            query_to_id_to_contexts[query] = id_to_contexts
        return query_to_id_to_contexts

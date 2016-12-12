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

def _create_kmer_index_dict(sequence_dictionary, k):
    """
    Return nested dictionary from kmer -> sequence ID -> position list
    """
    kmer_index_dict = defaultdict(lambda: defaultdict(list))
    for (key, sequence) in sequence_dictionary.items():
        for i in range(len(sequence) - k):
            kmer = sequence[i:i + k]
            kmer_index_dict[kmer][key].append(i)
    return kmer_index_dict


def _convert_nested_dict_values_to_to_sets(nested_dict):
    """
    Given a nested dictionary from kmer -> ID -> list of positions
    convert it to a dictionary from kmer -> ID -> set of positions.

    This simplifies code which wants to use these positions as a set
    against which membership is repeatedly checked.
    """
    return {
        kmer: {
            key: set(position_list)
            for (key, position_list) in id_to_positions_dict.items()}
        for (kmer, id_to_positions_dict)
        in nested_dict.items()
    }

class KmerLookup(object):
    def __init__(self, sequence_dictionary, min_kmer_size):
        """
        Simple initializer saves arguments but doesn't actually do any work,
        you must run the index() method.
        """
        self.sequence_dictionary = sequence_dictionary
        self.min_kmer_size = min_kmer_size
        self._kmer_index_dict = None

    def index(self):
        """
        Create index of kmer positions in the sequence dictionary supplied
        to the initializer.
        """
        self._kmer_index_dict = _create_kmer_index_dict(
            self.sequence_dictionary,
            self.min_kmer_size)
        self._kmer_index_dict = _convert_nested_dict_values_to_to_sets(
            self._kmer_index_dict)
        return self

    def _check_index_created(self):
        if self._kmer_index_dict is None:
            raise ValueError(
                "Index not created, you must call %s.index()" %
                self.__class__.__name__)
        if len(self._kmer_index_dict) == 0:
            raise ValueError("Kmer index is empty")

    def find_occurrences(self, query):
        """
        Returns a dictionary of sequence identifiers mapping to sets of
        offsets where the query subsequence occurs.
        """
        self._check_index_created()

        k = self.min_kmer_size
        kmer_index_dict = self._kmer_index_dict
        if len(query) < k:
            raise ValueError(
                "Query sequence '%s' must be at least %d characters" % (
                    query, k))

        first_kmer = query[:k]
        if first_kmer not in kmer_index_dict:
            return {}

        hits = kmer_index_dict[first_kmer]

        for i in range(1, len(query) - k):
            new_hits = defaultdict(set)
            kmer = query[i:i + k]
            for identifier, offset in kmer_index_dict.get(kmer, []):
                if identifier in hits and offset - 1 in hits[identifier]:
                    new_hits[identifier].add(offset)
            if len(new_hits) == 0:
                break
            hits = new_hits
        return {
            identifier: {offset - len(query) + k for offset in offsets}
            for (identifier, offsets) in hits.items()
        }

"""See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


Representation of an exon.

Examples:
    intron = Intron(1, 10, "+", "X")
    intron.getSequence()
"""

from sequence import Sequence


class Intron:

    canonical_splice_sites = ["GTAG", "ATAC", "GCAG"]
    fasta_file = None

    def __init__(
        self,
        exons: list,
        fasta_file: str = None,
        public_identifier: str = None,
    ) -> None:
        self.build_intron(exons)
        if fasta_file is not None:
            self.fasta_file = fasta_file
        self.sequence = None
        self.public_identifier = public_identifier

    def build_intron(self, exons: list) -> None:
        if exons[0].start > exons[1].start:
            exons.sort(key=lambda x: x.start)
            print("Left exon start coord > right exon start coord, will swap")

        exon_left = exons[0]
        exon_right = exons[1]

        self.start = exon_left.end + 1
        self.end = exon_right.start - 1
        self.length = self.end - self.start + 1
        self.strand = exon_left.strand
        self.location_name = exon_left.location_name
        if hasattr(exon_left, "fasta_file"):
            self.fasta_file = exon_left.fasta_file

    def get_sequence(self) -> str:
        if self.sequence is None:
            sequence = Sequence(
                self.start, self.end, self.strand, self.location_name, self.fasta_file
            )
            self.sequence = sequence

        if self.sequence.sequence is None:
            sequence.get_sequence()

        return self.sequence.sequence

    def is_splice_canonical(self) -> bool:
        sequence = self.get_sequence()
        donor = sequence[:2]
        acceptor = sequence[-2:]
        splice_site = donor + acceptor
        if splice_site in Intron.canonical_splice_sites:
            return True
        else:
            return False

    def intron_string(self, verbose: bool = False) -> str:
        start = self.start
        end = self.end
        if self.strand == "-":
            start = self.end
            end = self.start

        intron_string = "<" + str(start) + ".." + str(end) + ">"

        if verbose:
            intron_string = intron_string + ":" + self.strand + ":" + self.location_name

        return intron_string

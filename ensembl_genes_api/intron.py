# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
Representation of an intron.

Examples:
    exon1 = Exon(1, 10, "+", "X")
    exon2 = Exon(136, 150, "+", "X")
    intron = Intron([exon1, exon2])
    intron.getSequence()
"""

from sequence import Sequence


class Intron:
    """Representation of an intron.

    The location of an intron is always provided 5' -> 3' relative to the cDNA.
    Thus the start should be greater than the end if the transcript is on the
    reverse strand.

    Attributes:
      start: Start position of the intron.
      end: End position of the intron.
      strand: Strand of the intron, either + or -.
      location_name: Name of the region the intron is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the intron.
      public_identifier: Name of the intron.
      exons: The exons on each side of the introns
      canonical_splice_sites: list of canonical splice site, i.e, GTAG,...

    Raises:
      Exception: when exons size is different from 2
    """

    # pylint: disable=too-many-instance-attributes
    canonical_splice_sites = ["GTAG", "ATAC", "GCAG"]
    fasta_file = None

    def __init__(
        self,
        exons: list,
        fasta_file: str = None,
        public_identifier: str = None,
    ) -> None:
        """Initialise intron

        Args:
          exons: The exons on each side of the introns
          fasta_file: Path to a FASTA file containing the DNA of the region.
          public_identifier: Name of the intron.

        Raises:
          Exception: when exons size is different from 2
        """
        if len(exons) != 2:
            raise Exception(f"exons should only have 2 elements, not {len(exons)}")
        self.build_intron(exons)
        if fasta_file is not None:
            self.fasta_file = fasta_file
        self.sequence = None
        self.public_identifier = public_identifier

    def build_intron(self, exons: list) -> None:
        """Populate the intron attributes using the exons provided.

        Args:
          exons: Two exons delimiting the intron.
        """
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
        """The DNA sequence of the intron from the FASTA file attached.

        Returns:
          The DNA sequence.
        """
        if self.sequence is None:
            self.sequence = Sequence(
                self.start, self.end, self.strand, self.location_name, self.fasta_file
            )

        if self.sequence.sequence is None:
            self.sequence.get_sequence()

        return self.sequence.sequence

    def is_splice_canonical(self) -> bool:
        """Assert the splice junction is canonical

        The canonical splice junction are the most common two amino acids at the beginning
        and at the end of the intron. Non canonical splice junctions occur but usually
        means that either there is an error in the assembly or there is a misalignment.

        Returns:
          True if the two first amino acids and the two last amino acids concatenated matches
          one of:
            * GTAG
            * ATAC
            * CGAG
          False in any other cases
        """
        sequence = self.get_sequence()
        if sequence is None:
            return False
        donor = sequence[:2]
        acceptor = sequence[-2:]
        splice_site = donor + acceptor
        return bool(splice_site in Intron.canonical_splice_sites)

    def intron_string(self, verbose: bool = False) -> str:
        """Unique intron identifier based on its location.

        Args:
          verbose: add the strand and region name to the returned string.

        Returns:
          The unique identifier as (start..end) based on its direction:
          <1..10> if it is on the forward strand (+)
          <10..1> if it is on the reverse strand (-)

          When verbose is set the format is <1..10>:1:X
        """
        start = self.start
        end = self.end
        if self.strand == "-":
            start = self.end
            end = self.start

        intron_string = f"<{str(start)}..{str(end)}>"

        if verbose:
            intron_string = f"{intron_string}:{self.strand}:{self.location_name}"

        return intron_string

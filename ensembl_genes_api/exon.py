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
Representation of an exon.

Examples:
    exon = Exon(1, 10, "+", "X")
    exon.get_sequence()
"""


from sequence import Sequence


class Exon:  # pylint: disable=too-many-instance-attributes
    """Representation of an exon.

    The location of an exon is always provided 5' -> 3' on the forward strand.
    Thus the start should always be lower than the end.

    We use phase and end phase to provide information about the coding potential
    and where the next codon start at the start of the exon. The explanation from
    the Ensembl API:
    The Ensembl phase convention can be thought of as
    "the number of bases of the first codon which are
    on the previous exon".  It is therefore 0, 1 or 2
    (or -1 if the exon is non-coding).  In ascii art,
    with alternate codons represented by B<###> and
    B<+++>:

           Previous Exon   Intron   This Exon
        ...-------------            -------------...

        5'                    Phase                3'
        ...#+++###+++###          0 +++###+++###+...
        ...+++###+++###+          1 ++###+++###++...
        ...++###+++###++          2 +###+++###+++...

    Here is another explanation from Ewan:

    Phase means the place where the intron lands
    inside the codon - 0 between  codons, 1 between
    the 1st and second base, 2 between the second and
    3rd  base. Exons therefore have a start phase and
    a end phase, but introns have just one phase.

    Attributes:
      start: Start position of the exon.
      end: End position of the exon.
      strand: Strand of the exon, either + or -.
      location_name: Name of the region the exon is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the exon.
      public_identifier: Name of the exon.
      exon_start_phase: Phase of the exon, either -1, 0, 1 or 2.
      exon_end_phase: End phase of the exon, either -1, 0, 1 or 2.

    Raises:
      Exception: when start is greater than end
    """

    # TODO: exon_start_phase and exon_end_phase could have at least exon_ removed
    def __init__(  # pylint: disable=too-many-arguments
        self,
        start: int,
        end: int,
        strand: str,
        location_name: str,
        fasta_file: str = None,
        sequence: str = None,
        public_identifier: str = None,
        exon_start_phase: int = None,
        exon_end_phase: int = None,
    ) -> None:
        """Initialise Exon.

        Args:
          start: Start position of the exon.
          end: End position of the exon.
          strand: Strand of the exon, either + or -.
          location_name: Name of the region the exon is on.
          fasta_file: Path to a FASTA file containing the DNA of the region.
          sequence: DNA sequence of the exon.
          public_identifier: Name of the exon.
          exon_start_phase: Phase of the exon, either -1, 0, 1 or 2.
          exon_end_phase: End phase of the exon, either -1, 0, 1 or 2.

        Raises:
          Exception: when start is greater than end.
        """
        if start > end:
            raise Exception(f"Exon start is greater then exon end: {start} > {end}")

        self.start = start
        self.end = end
        self.strand = strand
        self.location_name = location_name
        self.fasta_file = fasta_file
        self.sequence = sequence
        self.public_identifier = public_identifier
        self.exon_start_phase = exon_start_phase
        self.exon_end_phase = exon_end_phase

    def get_sequence(self) -> str:
        """The DNA sequence of the exon from the FASTA file attached.

        Returns:
          The DNA sequence as a string.
        """

        if self.sequence is None:
            sequence = Sequence(
                self.start, self.end, self.strand, self.location_name, self.fasta_file
            )
            self.sequence = sequence

        if self.sequence.sequence is None:
            sequence.get_sequence()

        return self.sequence.sequence

    # TODO: this could be __str__ so it could simply be called in "String" context
    def exon_string(self, verbose: bool = False) -> str:
        """Unique exon identifier based on its location.

        Args:
          verbose: add the strand and region name to the returned string.

        Returns:
          The unique identifier as (start..end) based on its direction:
          (1..10) if it is on the forward strand (+)
          (10..1) if it is on the reverse strand (-)

          When verbose is set the format is (1..10):1:X
        """

        start = self.start
        end = self.end
        if self.strand == "-":
            start = self.end
            end = self.start

        exon_string = f"({start}..{end})"

        if verbose:
            exon_string = f"{exon_string}:{self.strand}:{self.location_name}"

        return exon_string

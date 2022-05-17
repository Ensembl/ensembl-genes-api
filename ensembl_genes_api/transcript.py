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
Representation of a transcript.

Examples:
    transcript = transcript(exons)
"""

import tempfile
import subprocess
import io
import os
import re
import copy

from intron import Intron


class Transcript:
    """Representation of a transcript.

    The location of a transcript is always provided 5' -> 3' on the forward strand.
    Thus the start should always be lower than the end.

    Attributes:
      translate_path: Path to the software 'translate' to translate DNA sequence into protein.
      exons: List of exons constituting the transcript.
      start: Start position of the transcript.
      end: End position of the transcript.
      strand: Strand of the transcript, either + or -.
      location_name: Name of the region the transcript is on.
      fasta_file: Path to a FASTA file containing the DNA of the region.
      sequence: DNA sequence of the transcript.
      public_identifier: Name of the transcript.
      introns: List of introns constituting the transcript, should be 0 if exons size is one.
      cds_genomic_start: Genomic start position of the CDS, always be lower than cds_genomic_end.
      cds_genomic_end: Genomic end position of the CDS, always be greater than cds_genomic_start.
      cds_sequence: Translateable cDNA sequence.
      translation_sequence: Protein sequence translated from cds_sequence.
    """

    # pylint: disable=too-many-instance-attributes
    translate_path = "translate"

    def __init__(
        self,
        exons: list,
        fasta_file: str = None,
        public_identifier: str = None,
    ) -> None:
        """Create the transcript object.

        Args:
          exons: List of exons which consitute the transcript.
          fasta_file: Path to a FASTA file containing the DNA of the region.
          public_identifier: Name of the transcript.
        """
        self.exons = exons
        self.fasta_file = fasta_file
        self.build_transcript(exons)
        self.sequence = None
        self.cds_genomic_start = None
        self.cds_genomic_end = None
        self.cds_sequence = None
        self.translation_sequence = None
        self.public_identifier = public_identifier

    @property
    def start(self) -> int:
        """Get :attr: start."""
        return self._start

    @property
    def end(self) -> int:
        """Get :attr: end."""
        return self._end

    @property
    def strand(self) -> str:
        """Get :attr: strand."""
        return self._strand

    @property
    def location_name(self) -> int:
        """Get :attr: location_name."""
        return self._location_name

    @property
    def introns(self) -> int:
        """Get :attr: introns."""
        return self._introns

    def build_transcript(self, exons: list) -> None:
        """Set the basic attributes of the transcript based on the exons.

        It will calculate the start and end of the transcript based on the first exon and last exon
        from exons.
        It will also generate the introns when the transcript has at least two exons.

        Args:
          exons: List of exons constituting the transcript.

        Raises:
          Exception: when exons have different strand.
          Exception: when exons have different location_name values.
          Exception: when start is greater than end
        """
        # Check the integrity of the exons
        strand = exons[0].strand
        location_name = exons[0].location_name
        for exon in exons:
            if exon.strand != strand:
                raise Exception(
                    "Inconsistent strands on the exons. Transplicing not supported"
                )
            if exon.location_name != location_name:
                raise Exception(
                    f"Exons should belong to the same sequence: "
                    f"{exon.location_name} {location_name}"
                )

        exons.sort(key=lambda x: x.start)
        self._start = exons[0].start
        self._end = exons[-1].end
        if strand == "-":
            exons.reverse()

        if self._start >= self._end:
            raise Exception("Transcript start was >= end, this should not be")

        # If we have multiple exons then calculate the introns
        self._introns = []
        if len(exons) > 1:
            for idx, exon in enumerate(exons[:-1]):
                intron = Intron([exon, exons[idx + 1]])
                self._introns.append(intron)

        self._strand = strand
        self._location_name = location_name

    def add_exons(self, exons: list) -> None:
        """Add exon(s) to the transcript.

        Args:
          exons: List of exon(s) to add to the transcript.
        """
        # Add a list of exons onto the existing set of exons. Rebuild transcript and
        # set any sequence to None, since that will need to be re-calculated
        self.exons = self.exons + exons
        self.build_transcript(self.exons)
        self.sequence = None

    def get_sequence(self) -> str:
        """Retrieve the DNA sequence of the transcript.

        Returns:
          The cDNA sequence.
        """
        if self.sequence is None:
            sequence = ""
            for exon in self.exons:
                sequence = sequence + exon.get_sequence()
            self.sequence = sequence

        return self.sequence

    def get_cds_sequence(self) -> str:
        """Retrieve the cDNA sequence of the translateable transcript sequence.

        Returns:
          The cDNA sequence without UTR(s).
        """
        if (
            self.cds_sequence is None
            and self.cds_genomic_start is not None
            and self.cds_genomic_end is not None
        ):
            self.construct_cds(self.cds_genomic_start, self.cds_genomic_end, self.exons)

        return self.cds_sequence

    def get_translation_sequence(self) -> str:
        """Retrieve the protein sequence of the transcript.

        Returns:
          The translation of the transcript.
        """
        cds_sequence = self.get_cds_sequence()
        if cds_sequence is not None:
            self.construct_translation(cds_sequence)

        return self.translation_sequence

    def construct_cds(self, genomic_start: int, genomic_end: int, exons: list) -> None:
        """Generate the CDS sequence.

        Args:
          genomic_start: Start of the CDS 5' -> 3' on the forward strand.
          genomic_end: End of the CDS 5' -> 3' on the forward strand.
          exons: List of exons.

        Raises:
          Exception: when start exon or end exon where not found
        """
        # Make a copy of the exons, based on start, then loop over the range of exons that
        # the cds start and end cover. Then edit the boundaries of the start and end exons
        # At the moment I've just made a temp transcript with the cds exons and directly
        # translates that
        forward_sorted_exons = copy.deepcopy(exons)
        forward_sorted_exons.sort(key=lambda x: x.start)
        start_exon_idx = Transcript.get_feature_index(
            genomic_start, forward_sorted_exons
        )
        end_exon_idx = Transcript.get_feature_index(genomic_end, forward_sorted_exons)

        if start_exon_idx is None or end_exon_idx is None:
            raise Exception(
                "Start or end exon index was not found based on genomic coordinate, this is wrong"
            )

        cds_exons = []
        for idx in range(start_exon_idx, end_exon_idx + 1):
            exon = forward_sorted_exons[idx]
            if idx == start_exon_idx and start_exon_idx == end_exon_idx:
                exon.sequence = None
                exon.start = genomic_start
                exon.end = genomic_end
                cds_exons.append(exon)
            elif idx == start_exon_idx:
                exon.sequence = None
                exon.start = genomic_start
                cds_exons.append(exon)
            elif idx == end_exon_idx:
                exon.sequence = None
                exon.end = genomic_end
                cds_exons.append(exon)
            else:
                cds_exons.append(exon)

        tmp_transcript = Transcript(cds_exons)
        cds_sequence = tmp_transcript.get_sequence().upper()

        self.cds_sequence = cds_sequence

    def construct_translation(self, cds_sequence: str) -> None:
        """Generate the protein sequence from the transcript.

        Args:
          cds_sequence: cDNA sequence of the translation.
        """
        # Just does a direct translation of a cds sequence that has already been calculated
        self.translation_sequence = Transcript.local_translate(cds_sequence)

    def compute_translation(self) -> None:
        """Generate the longest translation for the transcript."""
        # First remove any existing cds/translation info
        self.cds_sequence = None
        self.cds_genomic_start = None
        self.cds_genomic_end = None
        self.translation_sequence = None
        translations_met_required = []
        translations_met_optional = []
        translations_met_required = Transcript.run_translate(self.get_sequence(), True)
        translations_met_optional = Transcript.run_translate(self.get_sequence(), False)

        best_translation_met = None
        best_translation_no_met = None
        primary_translation = None

        if len(translations_met_required) > 0:
            best_translation_met = translations_met_required[0]

        if len(translations_met_optional) > 0:
            best_translation_no_met = translations_met_optional[0]

        if best_translation_met and best_translation_no_met:
            if translations_met_required[0][2] < 100 and translations_met_optional[0][
                2
            ] > (2 * translations_met_required[0][2]):
                primary_translation = best_translation_no_met
            else:
                primary_translation = best_translation_met
        elif best_translation_met:
            primary_translation = best_translation_met
        elif best_translation_no_met:
            primary_translation = best_translation_no_met

        if primary_translation is not None:
            sequence_start = primary_translation[0]
            sequence_end = primary_translation[1]
            if self.strand == "+":
                self.cds_genomic_start = Transcript.sequence_to_genomic_coord(
                    sequence_start, self.exons
                )
                self.cds_genomic_end = Transcript.sequence_to_genomic_coord(
                    sequence_end, self.exons
                )
            else:
                self.cds_genomic_start = Transcript.sequence_to_genomic_coord(
                    sequence_end, self.exons
                )
                self.cds_genomic_end = Transcript.sequence_to_genomic_coord(
                    sequence_start, self.exons
                )

            # Now store the cds and translation sequence
            self.get_translation_sequence()

    @staticmethod
    def run_translate(
        sequence, require_methionine: bool = False, min_length: int = 50
    ) -> list:
        """Run the translate software to find all possible ORFs.

        Args:
          require_methionine: Set to True if all possible ORFs need to start with a methionine.
          min_length: Minimum length of an ORF to be reported, set to 0 to disable.

        Returns:
          A list of a list of three elements:
          - start
          - end
          - length
        """
        translations = []
        translate_path = Transcript.translate_path
        with tempfile.NamedTemporaryFile(
            mode="w+t", delete=False
        ) as sequence_temp_file:
            sequence_temp_file.write(">tempseq\n" + sequence + "\n")
            sequence_temp_file.close()

            translate_command = [translate_path]
            if require_methionine:
                translate_command.append("-m")

            if min_length > 1:
                translate_command.append("-l")
                translate_command.append(str(min_length))

            translate_command.append(sequence_temp_file.name)

            translate_output = subprocess.Popen(
                translate_command, stdout=subprocess.PIPE
            )

            fasta_regex = re.compile(r"^>.+length (\d+),\s+nt (\d+)\.\.(\d+)")
            for line in io.TextIOWrapper(translate_output.stdout, encoding="utf-8"):
                match = fasta_regex.search(line)
                if match:
                    start = int(match.group(2))
                    end = int(match.group(3))
                    if start < end:
                        translations.append([start, end, int(match.group(1))])

            os.remove(sequence_temp_file.name)

        return translations

    @staticmethod
    def get_feature_index(genomic_position: int, features: list) -> int:
        """Return the index of the first feature which overlaps the genomic position.

        Args:
          genomic_position: A genomic position.
          features: List of features which may contain a feature overlapping the genomic position.

        Returns:
          The index of the first feature overlapping the genomic position or None.
        """
        for idx, feature in enumerate(features):
            if feature.start <= genomic_position <= feature.end:
                return idx
        return None

    @staticmethod
    def sequence_to_genomic_coord(
        sequence_position: int,
        features: list,
    ) -> int:
        """Convert the position from sequence coordinates to genomic coordinates.

        Args:
          sequence_position: Position to convert.

        Returns:
          The genomic position or None.
        """
        # This loops through a set features with an associated sequence to place a pair of sequence
        # coords onto the genome
        # A couple of straightforward use cases are converting protein and transcript coords to
        # genomic coords
        # Since the sequence might not cover all features, a feature start and end offset can be
        # provided
        # For example a CDS sequence might start/end in the middle of an exon and thus the offset
        # is needed
        combined_length = 0
        for feature in features:
            next_combined_length = combined_length + (feature.end - feature.start + 1)
            if next_combined_length >= sequence_position:
                remaining_offset = (sequence_position - combined_length) - 1
                if feature.strand == "+":
                    return feature.start + remaining_offset
                return feature.end - remaining_offset

            combined_length = next_combined_length

        return None

    @staticmethod
    def local_translate(sequence: str) -> str:
        """Translate the sequence in pure Python.

        Args:
          sequence: cDNA sequence to translate.

        Raises:
          Exception: when the length of sequence is not a multiple of 3.

        Returns:
          The protein sequence generated.
        """
        # fmt: off
        translation_table = {
            "AAA": "K", "CAA": "Q", "GAA": "E", "TAA": "*",
            "AAC": "N", "CAC": "H", "GAC": "D", "TAC": "Y",
            "AAG": "K", "CAG": "Q", "GAG": "E", "TAG": "*",
            "AAT": "N", "CAT": "H", "GAT": "D", "TAT": "Y",
            "ACA": "T", "CCA": "P", "GCA": "A", "TCA": "S",
            "ACC": "T", "CCC": "P", "GCC": "A", "TCC": "S",
            "ACG": "T", "CCG": "P", "GCG": "A", "TCG": "S",
            "ACT": "T", "CCT": "P", "GCT": "A", "TCT": "S",
            "AGA": "R", "CGA": "R", "GGA": "G", "TGA": "*",
            "AGC": "S", "CGC": "R", "GGC": "G", "TGC": "C",
            "AGG": "R", "CGG": "R", "GGG": "G", "TGG": "W",
            "AGT": "S", "CGT": "R", "GGT": "G", "TGT": "C",
            "ATA": "I", "CTA": "L", "GTA": "V", "TTA": "L",
            "ATC": "I", "CTC": "L", "GTC": "V", "TTC": "F",
            "ATG": "M", "CTG": "L", "GTG": "V", "TTG": "L",
            "ATT": "I", "CTT": "L", "GTT": "V", "TTT": "F",
        }
        # fmt: on

        translation = ""
        if len(sequence) % 3 == 0:
            for i in range(0, len(sequence), 3):
                codon = sequence[i : i + 3]
                translation += translation_table[codon]

        else:
            raise Exception(
                "Sequence passed in for local translation was not zero mod three"
            )

        return translation

    def transcript_string(self) -> str:
        """Unique string representation of Transcript.

        Generate a string representing the transcript uniquely by using the location of the
        transcript but also of the exon(s) and intron(s).

        Returns:
          The string representing the transcript in the following format:
          transcript; location='X'; strand='+'; structure=(1..30)<31..55>(56..102)
          transcript; location='X'; strand='-'; structure=(102..56)<55..31>(30..1)
        """
        transcript_string = (
            "transcript; location='"
            + self.location_name
            + "'; strand='"
            + self.strand
            + "'; structure="
        )
        intron_count = len(self.exons) - 1
        for idx, exon in enumerate(self.exons):
            transcript_string = transcript_string + exon.exon_string()
            if idx < intron_count:
                transcript_string = (
                    transcript_string + self.introns[idx].intron_string()
                )

        return transcript_string

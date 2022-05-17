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
Tests for the transcript module.
"""
import pytest

from exon import Exon
from transcript import Transcript

# This is needed because I use fixtures and I need to provide the function name as variable
# pylint: disable=redefined-outer-name.


@pytest.fixture
def example_transcript_forward() -> Transcript:
    """Create a simple Transcript on the forward strand for tests.

    Returns:
        Transcript object.
    """
    exon1 = Exon(1, 19, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa")
    exon2 = Exon(
        72, 105, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa"
    )
    return Transcript([exon1, exon2])


@pytest.fixture
def example_transcript_reverse() -> Transcript:
    """Create a simple Transcript on the reverse strand for tests.

    Returns:
        Transcript object.
    """
    exon1 = Exon(
        154, 170, "-", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa"
    )
    exon2 = Exon(19, 70, "-", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa")
    return Transcript([exon1, exon2])


@pytest.fixture
def example_transcript_sequence() -> Transcript:
    """Create a simple Transcript to test translations.

    Returns:
        Transcript object.
    """
    exon1 = Exon(1, 19, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa")
    exon2 = Exon(
        72, 106, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa"
    )
    exon3 = Exon(
        145, 185, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa"
    )
    exon4 = Exon(
        220, 260, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa"
    )
    exon5 = Exon(
        310, 370, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa"
    )
    return Transcript([exon1, exon2, exon3, exon4, exon5])


def test_get_start(example_transcript_forward: Transcript):
    """Test :attr: start getter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert example_transcript_forward.start == 1


def test_get_end(example_transcript_forward: Transcript):
    """Test :attr: end getter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert example_transcript_forward.end == 105


def test_get_strand(example_transcript_forward: Transcript):
    """Test :attr: strand getter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert example_transcript_forward.strand == "+"


def test_get_location(example_transcript_forward: Transcript):
    """Test :attr: location_name getter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert example_transcript_forward.location_name == "X"


def test_get_id(example_transcript_forward: Transcript):
    """Test :attr: public_identifier getter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert example_transcript_forward.public_identifier is None


def test_get_fasta_file(example_transcript_forward: Transcript):
    """Test :attr: fasta_file getter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert example_transcript_forward.fasta_file is None


def test_get_introns(example_transcript_forward: Transcript):
    """Test :attr: introns getter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert len(example_transcript_forward.introns) == 1


def test_get_introns_none():
    """Test :attr: introns getter when transcript has one exon."""
    exon = Exon(72, 105, "+", "X")
    transcript = Transcript([exon])
    assert len(transcript.introns) == 0


def test_set_id(example_transcript_forward: Transcript):
    """Test :attr: public_identifier setter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    example_transcript_forward.public_identifier = "transcript1"
    assert example_transcript_forward.public_identifier == "transcript1"


def test_set_fasta_file(example_transcript_forward: Transcript):
    """Test :attr: fasta_file setter.

    Args:
        example_transcript_forward: an Transcript object.
    """
    example_transcript_forward.fasta_file = "tests/data_genomic.fa"
    assert example_transcript_forward.fasta_file == "tests/data_genomic.fa"


def test_transcript_string_forward(example_transcript_forward: Transcript):
    """Test :func: transcript_string without verbose set on a forward object.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert (
        example_transcript_forward.transcript_string()
        == "transcript; location='X'; strand='+'; structure=(1..19)<20..71>(72..105)"
    )


def test_get_sequence(example_transcript_forward: Transcript):
    """Test :func: get_sequence() returns the DNA sequence on a forward object.

    Args:
        example_transcript_forward: an Transcript object.
    """
    assert (
        example_transcript_forward.get_sequence()
        == "ATGATATCGTATCATGATCTCCCAGCTCCAGTCAGTCAGCCAAGTCCCAGCTC"
    )


def test_add_exon(example_transcript_forward: Transcript):
    """Test :func: add_exon on a forward transcript.

    Args:
        example_transcript_forward: an Transcript object.
    """
    example_transcript_forward.add_exons([Exon(129, 135, "+", "X")])
    assert (
        example_transcript_forward.transcript_string()
        == "transcript; location='X'; strand='+'; "
        + "structure=(1..19)<20..71>(72..105)<106..128>(129..135)"
    )


def test_add_exon_reverse(example_transcript_reverse: Transcript):
    """Test :func: add_exon on a reverse transcript.

    Args:
        example_transcript_reverse: an Transcript object.
    """
    example_transcript_reverse.add_exons([Exon(1, 9, "-", "X")])
    assert (
        example_transcript_reverse.transcript_string()
        == "transcript; location='X'; strand='-'; "
        + "structure=(170..154)<153..71>(70..19)<18..10>(9..1)"
    )


def test_transcript_string_reverse(example_transcript_reverse: Transcript):
    """Test :func: transcript_string on a reverse object.

    Args:
        example_transcript_reverse: an Transcript object.
    """
    assert (
        example_transcript_reverse.transcript_string()
        == "transcript; location='X'; strand='-'; "
        + "structure=(170..154)<153..71>(70..19)"
    )


def test_get_sequence_reverse(example_transcript_reverse: Transcript):
    """Test :func: get_sequence() returns the DNA sequence on a reverse object.

    Args:
        example_transcript_reverse: an Transcript object.
    """
    assert (
        example_transcript_reverse.get_sequence()
        == "TCCTGTATAACGCTGTCTGAGTCGGATTAGCTATACGAACTGACATGATCATTGATATCGAACGACACG"
    )


def test_compute_translation_cds(example_transcript_sequence: Transcript):
    """Test :func: compute_translation() returns the longest protein sequence.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    example_transcript_sequence.compute_translation()
    assert (
        example_transcript_sequence.cds_sequence
        == "ATGATCTCCCAGCTCCAGTCAGTCAGCCAAGTCCCAGCTCACAGATATATGACAGCGTTA"
        + "TACAGGATAGGATCCAAAGTATTCATGTCAGTTCGTATCGCTAATCCGACTCAGTCCCAG"
        + "CTCTATCAGTATCGTGTCGTTCGATATCAAGTATCATGTCAGTTCGTATTGCTAAGTATA"
        + "TCG"
    )


def test_cds_genomic_start(example_transcript_sequence: Transcript):
    """Test :attr: cds_genomic_start getter.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    example_transcript_sequence.compute_translation()
    assert example_transcript_sequence.cds_genomic_start == 14


def test_cds_genomic_end(example_transcript_sequence: Transcript):
    """Test :attr: cds_genomic_end getter.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    example_transcript_sequence.compute_translation()
    assert example_transcript_sequence.cds_genomic_end == 369


def test_get_translation_sequence(example_transcript_sequence: Transcript):
    """Test :func: get_translation_sequence() returns the longest protein sequence.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    example_transcript_sequence.compute_translation()
    assert (
        example_transcript_sequence.get_translation_sequence()
        == "MISQLQSVSQVPAHRYMTALYRIGSKVFMSVRIANPTQSQLYQYRVVRYQVSCQFVLLSIS"
    )


def test_local_translate(example_transcript_sequence: Transcript):
    """Test :func: local_translate() returns a protein sequence.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    example_transcript_sequence.compute_translation()
    assert (
        Transcript.local_translate(example_transcript_sequence.get_cds_sequence())
        == "MISQLQSVSQVPAHRYMTALYRIGSKVFMSVRIANPTQSQLYQYRVVRYQVSCQFVLLSIS"
    )


def test_run_translate(example_transcript_sequence: Transcript):
    """Test :func: run_translate() returns the longest possible translated sequence.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    sequence = Transcript.run_translate(example_transcript_sequence.get_sequence())
    assert sequence[0] == [5, 196, 64]


def test_run_translate_meth(example_transcript_sequence: Transcript):
    """Test :func: run_translate() returns the translated sequence starting with M.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    sequence = Transcript.run_translate(
        example_transcript_sequence.get_sequence(), require_methionine=True
    )
    assert sequence[0] == [14, 196, 61]


def test_sequence_to_genomic_coord(example_transcript_reverse: Transcript):
    """Test :func: sequence_to_genomic_coord() returns a cDNA position in genomic coord.

    Args:
        example_transcript_reverse: an Transcript object.
    """
    assert (
        Transcript.sequence_to_genomic_coord(14, example_transcript_reverse.exons)
        == 157
    )


def test_get_feature_index(example_transcript_sequence: Transcript):
    """Test :func: get_feature_index() returns a cDNA position in genomic coord.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    assert Transcript.get_feature_index(149, example_transcript_sequence.exons) == 2


def test_get_feature_index_notfound(example_transcript_sequence: Transcript):
    """Test :func: get_feature_index() returns a cDNA position in genomic coord.

    Args:
        example_transcript_sequence: an Transcript object.
    """
    assert Transcript.get_feature_index(143, example_transcript_sequence.exons) is None


def test_fail_set_start(example_transcript_forward: Transcript):
    """Test :attr: start setter fails.

    Args:
        example_transcript_forward: an Transcript object.
    """
    with pytest.raises(AttributeError, match=r"can't set attribute"):
        example_transcript_forward.start = 2


def test_fail_set_end(example_transcript_forward: Transcript):
    """Test :attr: end setter fails.

    Args:
        example_transcript_forward: an Transcript object.
    """
    with pytest.raises(AttributeError, match=r"can't set attribute"):
        example_transcript_forward.end = 28


def test_fail_set_strand(example_transcript_forward: Transcript):
    """Test :attr: strand setter fails.

    Args:
        example_transcript_forward: an Transcript object.
    """
    with pytest.raises(AttributeError, match=r"can't set attribute"):
        example_transcript_forward.strand = "-"


def test_fail_set_location(example_transcript_forward: Transcript):
    """Test :attr: location_name setter fails.

    Args:
        example_transcript_forward: an Transcript object.
    """
    with pytest.raises(AttributeError, match=r"can't set attribute"):
        example_transcript_forward.location_name = "2"


def test_fail_different_sequence():
    """Test creating object when one of the exon is on a different sequence."""
    with pytest.raises(
        Exception, match=r"Exons should belong to the same sequence: 2 X"
    ):
        Transcript(
            [Exon(3, 20, "+", "X"), Exon(30, 42, "+", "2"), Exon(50, 62, "+", "X")]
        )


def test_fail_transsplicing():
    """Test creating object when one of the exon is on a different strand."""
    with pytest.raises(
        Exception,
        match=r"Inconsistent strands on the exons. Transplicing not supported",
    ):
        Transcript(
            [Exon(3, 20, "+", "X"), Exon(30, 42, "-", "X"), Exon(50, 62, "+", "X")]
        )


def test_fail_local_translate(example_transcript_forward: Transcript):
    """Test :func: local_translate() fails if sequence is not a a multiple of 3.

    Args:
        example_transcript_forward: an Transcript object.
    """
    with pytest.raises(
        Exception,
        match=r"Sequence passed in for local translation was not zero mod three",
    ):
        Transcript.local_translate(example_transcript_forward.get_sequence())

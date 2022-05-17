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
Tests for the intron module.
"""


import pytest

from exon import Exon
from intron import Intron

# This is needed because I use fixtures and I need to provide the function name as variable
# pylint: disable=redefined-outer-name


@pytest.fixture
def example_intron_forward() -> Intron:
    """Create a simple Intron on the forward strand for tests

    Returns:
        Intron object
    """
    exon1 = Exon(1, 19, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa")
    exon2 = Exon(
        72, 105, "+", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa"
    )
    return Intron([exon1, exon2])


@pytest.fixture
def example_intron_reverse() -> Intron:
    """Create a simple Intron on the reverse strand for tests

    Returns:
        Intron object
    """
    exon1 = Exon(
        154, 170, "-", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa"
    )
    exon2 = Exon(19, 70, "-", "X", fasta_file="ensembl_genes_api/tests/data_genomic.fa")
    return Intron([exon1, exon2])


def test_get_start(example_intron_forward: Intron):
    """Test :attr: start getter

    Args:
        example_intron_forward: an Intron object
    """
    assert example_intron_forward.start == 20


def test_get_end(example_intron_forward: Intron):
    """Test :attr: end getter

    Args:
        example_intron_forward: an Intron object
    """
    assert example_intron_forward.end == 71


def test_get_strand(example_intron_forward: Intron):
    """Test :attr: strand getter

    Args:
        example_intron_forward: an Intron object
    """
    assert example_intron_forward.strand == "+"


def test_get_location(example_intron_forward: Intron):
    """Test :attr: location_name getter

    Args:
        example_intron_forward: an Intron object
    """
    assert example_intron_forward.location_name == "X"


def test_get_id(example_intron_forward: Intron):
    """Test :attr: public_identifier getter

    Args:
        example_intron_forward: an Intron object
    """
    assert example_intron_forward.public_identifier is None


def test_get_fasta_file(example_intron_forward: Intron):

    """Test :attr: fasta_file getter

    Args:
        example_intron_forward: an Exon object
    """
    assert (
        example_intron_forward.fasta_file == "ensembl_genes_api/tests/data_genomic.fa"
    )


def test_set_start(example_intron_forward: Intron):
    """Test :attr: start setter

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.start = 2
    assert example_intron_forward.start == 2


def test_set_end(example_intron_forward: Intron):
    """Test :attr: end setter

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.end = 28
    assert example_intron_forward.end == 28


def test_set_strand(example_intron_forward: Intron):
    """Test :attr: strand setter

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.strand = "-"
    assert example_intron_forward.strand == "-"


def test_set_location(example_intron_forward: Intron):
    """Test :attr: location_name setter

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.location_name = "2"
    assert example_intron_forward.location_name == "2"


def test_set_id(example_intron_forward: Intron):
    """Test :attr: public_identifier setter

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.public_identifier = "intron1"
    assert example_intron_forward.public_identifier == "intron1"


def test_set_fasta_file(example_intron_forward: Intron):
    """Test :attr: fasta_file setter

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.fasta_file = "tests/data_genomic.fa"
    assert example_intron_forward.fasta_file == "tests/data_genomic.fa"


def test_intron_string_forward(example_intron_forward: Intron):
    """Test :func: intron_string without verbose set on a forward object

    Args:
        example_intron_forward: an Intron object
    """
    assert example_intron_forward.intron_string() == "<20..71>"


def test_intron_string_forward_verbose(
    example_intron_forward: Intron,
):  # pylint: disable=invalid-name
    """Test :func: intron_string with verbose set on a forward object

    Args:
        example_intron_forward: an Intron object
    """
    assert example_intron_forward.intron_string(True) == "<20..71>:+:X"


def test_get_sequence(example_intron_forward: Intron):
    """Test :func: get_sequence() returns the DNA sequence on a forward object

    Args:
        example_intron_forward: an Intron object
    """
    assert (
        example_intron_forward.get_sequence()
        == "GTGTCGTTCGATATCAATGATCATGTCAGTTCGTATAGCTAATCCGACTCAG"
    )


def test_is_splice_canonical(example_intron_forward: Intron):
    """Test :func: is_splice_canonical() returns True on a canonical splice site

    Args:
        example_intron_forward: an Intron object
    """
    assert example_intron_forward.is_splice_canonical()


def test_intron_string_forward_mod(example_intron_forward: Intron):
    """Test :func: intron_string without verbose set on a modified forward object

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.start = 3
    assert example_intron_forward.intron_string() == "<3..71>"


def test_get_sequence_mod(example_intron_forward: Intron):
    """Test :func: get_sequence() returns the DNA sequence on a modified forward object

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.start = 4
    assert (
        example_intron_forward.get_sequence()
        == "ATATCGTATCATGATCGTGTCGTTCGATATCAATGATCATGTCAGTTCGTATAGCTAATCCGACTCAG"
    )


def test_is_splice_canonical_mod(example_intron_forward: Intron):
    """Test :func: is_splice_canonical() returns False on a non canonical splice site

    Args:
        example_intron_forward: an Intron object
    """
    example_intron_forward.start = 4
    assert not example_intron_forward.is_splice_canonical()


def test_intron_string_reverse(example_intron_reverse: Intron):
    """Test :func: intron_string with verbose set on a reverse object

    Args:
        example_intron_reverse: an Intron object
    """
    assert example_intron_reverse.intron_string() == "<153..71>"


def test_intron_string_reverse_verbose(
    example_intron_reverse: Intron,
):  # pylint: disable=invalid-name
    """Test :func: intron_string with verbose set on a reverse object

    Args:
        example_intron_reverse: an Intron object
    """
    assert example_intron_reverse.intron_string(True) == "<153..71>:-:X"


def test_get_sequence_reverse(example_intron_reverse: Intron):
    """Test :func: get_sequence() returns the DNA sequence on a reverse object

    Args:
        example_intron_reverse: an Intron object
    """
    assert (
        example_intron_reverse.get_sequence()
        == "ATATATCTGTGTCATATCTGTTACATCGCGTATCTCTCTGCTCTGGCTGAGCTGGGACTTGGCTGACTGACTGGAGCTGGGAC"
    )


def test_is_splice_canonical_reverse(
    example_intron_reverse: Intron,
):  # pylint: disable=invalid-name
    """Test :func: is_splice_canonical() returns True on a canonical splice site from reverse

    Args:
        example_intron_reverse: an Intron object
    """
    assert example_intron_reverse.is_splice_canonical()


def test_fail_one_exon():
    """Test creating object when exons size is 1"""
    with pytest.raises(Exception, match=r"exons should only have 2 elements, not 1"):
        Intron([Exon(3, 20, "+", "X")])


def test_fail_three_exon():
    """Test creating object when exons size is 3"""
    with pytest.raises(Exception, match=r"exons should only have 2 elements, not 3"):
        Intron([Exon(3, 20, "+", "X"), Exon(30, 42, "+", "X"), Exon(50, 62, "+", "X")])

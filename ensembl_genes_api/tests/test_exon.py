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
limitations under the License."""
import pytest

from exon import Exon

# This is needed because I use fixtures and I need to provide the function name as variable
# pylint: disable=redefined-outer-name


@pytest.fixture
def example_exon_forward() -> Exon:
    """Create a simple Exon on the forward strand for tests

    Returns:
        Exon object"""
    return Exon(1, 30, "+", "X", fasta_file="ensembl_genes_api/tests/data_exon.fa")


@pytest.fixture
def example_exon_reverse() -> Exon:
    """Create a simple Exon on the reverse strand for tests

    Returns:
        Exon object"""
    return Exon(1, 30, "-", "X", fasta_file="ensembl_genes_api/tests/data_exon.fa")


def test_get_start(example_exon_forward: Exon):
    """Test :attr: start getter

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.start == 1


def test_get_end(example_exon_forward: Exon):
    """Test :attr: end getter

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.end == 30


def test_get_strand(example_exon_forward: Exon):
    """Test :attr: strand getter

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.strand == "+"


def test_get_location(example_exon_forward: Exon):
    """Test :attr: location_name getter

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.location_name == "X"


def test_get_start_phase(example_exon_forward: Exon):
    """Test :attr: exon_start_phase getter

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.exon_start_phase == -1


def test_get_end_phase(example_exon_forward: Exon):
    """Test :attr: exon_end_phase getter

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.exon_end_phase == -1


def test_get_id(example_exon_forward: Exon):
    """Test :attr: public_identifier getter

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.public_identifier is None


def test_get_fasta_file(example_exon_forward: Exon):

    """Test :attr: fasta_file getter

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.fasta_file == "ensembl_genes_api/tests/data_exon.fa"


def test_set_start(example_exon_forward: Exon):
    """Test :attr: start setter

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.start = 2
    assert example_exon_forward.start == 2


def test_set_end(example_exon_forward: Exon):
    """Test :attr: end setter

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.end = 28
    assert example_exon_forward.end == 28


def test_set_strand(example_exon_forward: Exon):
    """Test :attr: strand setter

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.strand = "-"
    assert example_exon_forward.strand == "-"


def test_set_location(example_exon_forward: Exon):
    """Test :attr: location_name setter

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.location_name = "2"
    assert example_exon_forward.location_name == "2"


def test_set_start_phase(example_exon_forward: Exon):
    """Test :attr: exon_start_phase setter

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.exon_start_phase = 0
    assert example_exon_forward.exon_start_phase == 0


def test_set_end_phase(example_exon_forward: Exon):
    """Test :attr: exon_end_phase setter

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.exon_end_phase = 1
    assert example_exon_forward.exon_end_phase == 1


def test_set_id(example_exon_forward: Exon):
    """Test :attr: public_identifier setter

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.public_identifier = "exon1"
    assert example_exon_forward.public_identifier == "exon1"


def test_set_fasta_file(example_exon_forward: Exon):
    """Test :attr: fasta_file setter

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.fasta_file = "tests/data_exon.fa"
    assert example_exon_forward.fasta_file == "tests/data_exon.fa"


def test_exon_string_forward(example_exon_forward: Exon):
    """Test :func: exon_string without verbose set on a forward object

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.exon_string() == "(1..30)"


def test_exon_string_forward_verbose(
    example_exon_forward: Exon,
):  # pylint: disable=invalid-name
    """Test :func: exon_string with verbose set on a forward object

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.exon_string(True) == "(1..30):+:X"


def test_get_sequence(example_exon_forward: Exon):
    """Test :func: get_sequence() returns the DNA sequence on a forward object

    Args:
        example_exon_forward: an Exon object
    """
    assert example_exon_forward.get_sequence() == "ATGATATCGTATCATGATCGTBTCGTTCGA"


def test_exon_string_forward_mod(example_exon_forward: Exon):
    """Test :func: exon_string without verbose set on a modified forward object

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.start = 3
    assert example_exon_forward.exon_string() == "(3..30)"


def test_get_sequence_mod(example_exon_forward: Exon):
    """Test :func: get_sequence() returns the DNA sequence on a modified forward object

    Args:
        example_exon_forward: an Exon object
    """
    example_exon_forward.start = 4
    assert example_exon_forward.get_sequence() == "ATATCGTATCATGATCGTBTCGTTCGA"


def test_exon_string_reverse(example_exon_reverse: Exon):
    """Test :func: exon_string with verbose set on a reverse object

    Args:
        example_exon_reverse: an Exon object
    """
    assert example_exon_reverse.exon_string() == "(30..1)"


def test_exon_string_reverse_verbose(
    example_exon_reverse: Exon,
):  # pylint: disable=invalid-name
    """Test :func: exon_string with verbose set on a reverse object

    Args:
        example_exon_reverse: an Exon object
    """
    assert example_exon_reverse.exon_string(True) == "(30..1):-:X"


def test_get_sequence_reverse(example_exon_reverse: Exon):
    """Test :func: get_sequence() returns the DNA sequence on a reverse object

    Args:
        example_exon_reverse: an Exon object
    """
    assert example_exon_reverse.get_sequence() == "TCGAACGABACGATCATGATACGATATCAT"


def test_fail_start_gt_end():
    """Test creating object when start is greater than end

    Args:
        example_exon_forward: an Exon object
    """
    with pytest.raises(Exception, match=r"Exon start is greater then exon end:"):
        Exon(30, 2, "+", "X")

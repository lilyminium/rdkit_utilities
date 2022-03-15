"""
Unit and regression test for the rdkit_utilities package.
"""

# Import package, test suite, and other packages as needed
import rdkit_utilities
import pytest
import sys

def test_rdkit_utilities_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "rdkit_utilities" in sys.modules

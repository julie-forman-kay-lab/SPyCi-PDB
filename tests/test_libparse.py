"""Test libparse."""
from spycipdb.libs import libparse


def test_values_to_dict():
    """Test list of values to dictionary."""
    input = ["par1=1", "par2='my name'"]
    expected = {'par1': 1, 'par2': 'my name'}
    
    test = libparse.values_to_dict(input)
    assert test == expected

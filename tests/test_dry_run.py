# tests/test_dry_run.py
import pytest
from bscampp.pipeline import bscampp_pipeline

def test_bscampp_pipeline():
    res = bscampp_pipeline(dry_run=True)
    assert res == True

# tests/test_dry_run.py
import pytest, os
from bscampp.pipeline import bscampp_pipeline

def test_bscampp_pipeline():
    res = bscampp_pipeline(dry_run=True)
    assert res == True

    # remove bscampp_output that's created
    if os.path.isdir('bscampp_output'):
        os.rmdir('bscampp_output')

# tests/test_dry_run.py
import pytest, os
from bscampp.pipeline import bscampp_pipeline, scampp_pipeline

# test BSCAMPP
def test_bscampp_pipeline():
    res = bscampp_pipeline(dry_run=True)
    assert res == True

    # remove bscampp_output that's created
    if os.path.isdir('bscampp_output'):
        os.rmdir('bscampp_output')

# test SCAMPP (almost the same as BSCAMPP)
def test_scampp_pipeline():
    res = scampp_pipeline(dry_run=True)
    assert res == True

    # remove scampp_output that's created
    if os.path.isdir('scampp_output'):
        os.rmdir('scampp_output')

"""Field types of the spatial viewer's bins payload (GET /api/runs/{rid}/bins): a positive_* column is a yes/no call
only when its values are yes/no; the share of depth-matching draws in which a bin is positive is a number."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from webapp import server


def _table() -> pd.DataFrame:
    n = 6
    return pd.DataFrame({
        "barcode": [f"b{k}" for k in range(n)], "sample_id": ["S1"] * n, "x_um": np.arange(n) * 16.0,
        "y_um": np.zeros(n), "tissue_class": ["stroma"] * n, "total_umi": np.arange(n) + 5,
        "in_region_of_interest": [True, False, True, False, True, False],
        "positive_osteoclast": [True, False, False, True, False, False],                   # booleans
        "positive_osteoblast": [1, 0, 0, 1, np.nan, 0],                                     # numbers 0 / 1
        "positive_granulocyte": ["True", "False", None, "False", "True", "False"],          # text (a CSV with gaps)
        "positive_fraction_osteoclast_pooled_depthmatched_draws": [0.0, 0.05, 1.0, 0.5, 0.0, 0.95],
    })


def test_is_yes_no():
    t = _table()
    assert server.is_yes_no(t["positive_osteoclast"])
    assert server.is_yes_no(t["positive_osteoblast"])
    assert server.is_yes_no(t["positive_granulocyte"])
    assert not server.is_yes_no(t["positive_fraction_osteoclast_pooled_depthmatched_draws"])
    assert server.is_yes_no(pd.Series([0.0, 1.0, 1.0]))   # a share that is 0 or 1 everywhere is a yes/no


def test_share_of_draws_is_a_number_and_calls_stay_yes_no():
    table = server.BinsTable(_table(), {})
    kinds = table.kinds
    assert kinds["positive_fraction_osteoclast_pooled_depthmatched_draws"] == "number"
    for col in ("in_region_of_interest", "positive_osteoclast", "positive_osteoblast", "positive_granulocyte"):
        assert kinds[col] == "bool", col
    fields = json.loads(table.payload("S1"))["fields"]
    assert fields["positive_fraction_osteoclast_pooled_depthmatched_draws"] == {
        "type": "number", "values": [0.0, 0.05, 1.0, 0.5, 0.0, 0.95]}
    assert fields["positive_osteoblast"] == {"type": "bool", "values": [1, 0, 0, 1, 0, 0]}
    assert fields["positive_granulocyte"]["values"] == [1, 0, 0, 0, 1, 0]
    assert fields["tissue_class"]["type"] == "category" and fields["total_umi"]["type"] == "number"


def test_share_of_one_draw_is_still_a_number():
    """With uncoupling.n_draws = 1 every share is 0.0 or 1.0; the column stays a number (Color by -> Numbers)."""
    t = _table()
    t["positive_fraction_osteoclast_pooled_depthmatched_draws"] = [0.0, 1.0, 1.0, 0.0, 0.0, 1.0]
    assert server.BinsTable(t, {}).kinds["positive_fraction_osteoclast_pooled_depthmatched_draws"] == "number"

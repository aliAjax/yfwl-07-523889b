from __future__ import annotations

import importlib


def pytest_runtest_logreport(report):
    if report.when != "call":
        return
    test_name = report.nodeid.rsplit("::", 1)[-1]
    mapping = getattr(importlib.import_module("test_outputs"), "TEST_CASE_IDS", {})
    outcome = "PASS" if report.passed else "FAIL"
    raw_ids = mapping.get(test_name, ())
    case_ids = (raw_ids,) if isinstance(raw_ids, str) else raw_ids
    for case_id in case_ids:
        print(f"CASE {case_id} {outcome}")

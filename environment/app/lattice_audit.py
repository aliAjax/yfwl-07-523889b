#!/usr/bin/env python3
"""Legacy modular residual auditor."""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any


LIMITS = {
    "min_observations": 6,
    "max_observations": 6,
    "max_dimension": 2,
    "min_moduli": 2,
    "max_moduli": 2,
    "max_coefficient": 2,
    "max_basis_entry": 8,
    "max_determinant": 32,
}


class InputProblem(ValueError):
    """The supplied modular ticket cannot be parsed."""


class DuplicateKey(InputProblem):
    """A JSON object repeats a key."""


def pairs_hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DuplicateKey(key)
        result[key] = value
    return result


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_bytes().decode("utf-8"),
            object_pairs_hook=pairs_hook,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise InputProblem("read") from exc
    if not isinstance(value, dict):
        raise InputProblem("object")
    return value


def exact_keys(value: dict[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise InputProblem("keys")


def whole(value: Any, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InputProblem("integer")
    if value < minimum or value > maximum:
        raise InputProblem("range")
    return value


def label(value: Any) -> str:
    if not isinstance(value, str) or not value:
        raise InputProblem("label")
    return value


def prime_test(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def validate(payload: dict[str, Any]) -> dict[str, Any]:
    exact_keys(
        payload,
        {"schema_version", "case_id", "coefficients", "moduli", "limits"},
    )
    whole(payload["schema_version"], 1, 1)
    case_id = label(payload["case_id"])

    limits = payload["limits"]
    if not isinstance(limits, dict):
        raise InputProblem("limits")
    exact_keys(limits, set(LIMITS))
    for key, expected in LIMITS.items():
        whole(limits[key], expected, expected)

    raw_coefficients = payload["coefficients"]
    if not isinstance(raw_coefficients, list):
        raise InputProblem("coefficients")
    if len(raw_coefficients) != LIMITS["max_observations"]:
        raise InputProblem("observation_count")
    coefficients = []
    for item in raw_coefficients:
        if not isinstance(item, list) or len(item) != 2:
            raise InputProblem("coefficient")
        first = whole(
            item[0],
            -LIMITS["max_coefficient"],
            LIMITS["max_coefficient"],
        )
        second = whole(
            item[1],
            -LIMITS["max_coefficient"],
            LIMITS["max_coefficient"],
        )
        if first == 0 and second == 0:
            raise InputProblem("zero_coefficient")
        coefficients.append((first, second))

    raw_moduli = payload["moduli"]
    if not isinstance(raw_moduli, list):
        raise InputProblem("moduli")
    if len(raw_moduli) != LIMITS["max_moduli"]:
        raise InputProblem("modulus_count")
    moduli = []
    seen: set[int] = set()
    for item in raw_moduli:
        if not isinstance(item, dict):
            raise InputProblem("modulus_object")
        exact_keys(item, {"modulus", "residuals"})
        modulus = whole(item["modulus"], 2, 2147483647)
        if not prime_test(modulus):
            raise InputProblem("modulus")
        if modulus in seen:
            raise InputProblem("duplicate_modulus")
        seen.add(modulus)
        raw_residuals = item["residuals"]
        if not isinstance(raw_residuals, list) or len(raw_residuals) != len(coefficients):
            raise InputProblem("residual_count")
        residuals = []
        for pair in raw_residuals:
            if not isinstance(pair, list) or len(pair) != 2:
                raise InputProblem("residual_pair")
            residuals.append(
                (
                    whole(pair[0], 0, modulus - 1),
                    whole(pair[1], 0, modulus - 1),
                )
            )
        moduli.append({"modulus": modulus, "residuals": residuals})
    return {
        "case_id": case_id,
        "coefficients": coefficients,
        "moduli": moduli,
    }


def build(parsed: dict[str, Any]) -> dict[str, Any]:
    observation_count = len(parsed["coefficients"])
    modulus_count = len(parsed["moduli"])
    return {
        "schema_version": 1,
        "case_id": parsed["case_id"],
        "basis": [[1, 0], [1, 1]],
        "check": {
            "observation_count": observation_count,
            "modulus_count": modulus_count,
            "satisfied_equations": 2 * observation_count * modulus_count,
        },
    }


def compact(value: dict[str, Any]) -> bytes:
    text_value = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return (text_value + "\n").encode("utf-8")


def install(path: Path, data: bytes, mode: int = 0o644) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary_path, mode)
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def execute(args: argparse.Namespace) -> None:
    parsed = validate(load_json(Path(args.input)))
    install(Path(args.output), compact(build(parsed)))


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args(argv)


def main() -> int:
    try:
        execute(parse_args(sys.argv[1:]))
    except Exception:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

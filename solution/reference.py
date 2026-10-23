#!/usr/bin/env python3
"""Recover a small primitive lattice basis from ordered modular equations."""

from __future__ import annotations

import argparse
import json
import math
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


class ContractError(ValueError):
    """The input does not match the public contract."""


class DuplicateKeyError(ContractError):
    """A JSON object repeats a key."""


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DuplicateKeyError(key)
        result[key] = value
    return result


def load_payload(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(
            path.read_bytes().decode("utf-8"),
            object_pairs_hook=_pairs,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ContractError("input") from exc
    if not isinstance(value, dict):
        raise ContractError("object")
    return value


def exact_keys(value: dict[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ContractError("keys")


def integer(value: Any, low: int, high: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractError("integer")
    if value < low or value > high:
        raise ContractError("range")
    return value


def text(value: Any) -> str:
    if not isinstance(value, str) or not value:
        raise ContractError("text")
    return value


def is_prime(value: int) -> bool:
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


def validate_limits(raw: Any) -> None:
    if not isinstance(raw, dict):
        raise ContractError("limits")
    exact_keys(raw, set(LIMITS))
    for key, expected in LIMITS.items():
        integer(raw[key], expected, expected)


def validate_coefficients(raw: Any) -> list[tuple[int, int]]:
    if not isinstance(raw, list):
        raise ContractError("coefficients")
    if len(raw) != LIMITS["min_observations"]:
        raise ContractError("observation_count")
    coefficients: list[tuple[int, int]] = []
    for item in raw:
        if not isinstance(item, list) or len(item) != 2:
            raise ContractError("coefficient")
        first = integer(item[0], -2, 2)
        second = integer(item[1], -2, 2)
        if first == 0 and second == 0:
            raise ContractError("zero_coefficient")
        coefficients.append((first, second))
    return coefficients


def validate_moduli(
    raw: Any, observation_count: int
) -> list[dict[str, Any]]:
    if not isinstance(raw, list):
        raise ContractError("moduli")
    if len(raw) != LIMITS["min_moduli"]:
        raise ContractError("modulus_count")
    moduli: list[dict[str, Any]] = []
    seen: set[int] = set()
    for item in raw:
        if not isinstance(item, dict):
            raise ContractError("modulus_object")
        exact_keys(item, {"modulus", "residuals"})
        modulus = integer(item["modulus"], 2, 2147483647)
        if not is_prime(modulus):
            raise ContractError("modulus")
        if modulus in seen:
            raise ContractError("duplicate_modulus")
        seen.add(modulus)
        raw_residuals = item["residuals"]
        if not isinstance(raw_residuals, list) or len(raw_residuals) != observation_count:
            raise ContractError("residual_count")
        residuals = []
        for pair in raw_residuals:
            if not isinstance(pair, list) or len(pair) != 2:
                raise ContractError("residual_pair")
            residuals.append(
                (
                    integer(pair[0], 0, modulus - 1),
                    integer(pair[1], 0, modulus - 1),
                )
            )
        moduli.append({"modulus": modulus, "residuals": residuals})
    return moduli


def validate(payload: dict[str, Any]) -> dict[str, Any]:
    exact_keys(
        payload,
        {"schema_version", "case_id", "coefficients", "moduli", "limits"},
    )
    integer(payload["schema_version"], 1, 1)
    case_id = text(payload["case_id"])
    validate_limits(payload["limits"])
    coefficients = validate_coefficients(payload["coefficients"])
    moduli = validate_moduli(payload["moduli"], len(coefficients))
    return {"case_id": case_id, "coefficients": coefficients, "moduli": moduli}


def solve_modular_system(
    coefficients: list[tuple[int, int]],
    residuals: list[int],
    modulus: int,
) -> tuple[int, int]:
    for left_index in range(len(coefficients)):
        first_a, first_b = coefficients[left_index]
        for right_index in range(left_index + 1, len(coefficients)):
            second_a, second_b = coefficients[right_index]
            determinant = (first_a * second_b - first_b * second_a) % modulus
            if determinant == 0:
                continue
            inverse = pow(determinant, -1, modulus)
            first_target = residuals[left_index]
            second_target = residuals[right_index]
            first = (
                (first_target * second_b - first_b * second_target) * inverse
            ) % modulus
            second = (
                (first_a * second_target - first_target * second_a) * inverse
            ) % modulus
            if all(
                (a_value * first + b_value * second - target) % modulus == 0
                for (a_value, b_value), target in zip(coefficients, residuals)
            ):
                return first, second
    raise ContractError("inconsistent_system")


def crt(residues: list[int], moduli: list[int]) -> int:
    product = math.prod(moduli)
    value = 0
    for residue, modulus in zip(residues, moduli):
        partial = product // modulus
        value += residue * partial * pow(partial, -1, modulus)
    return value % product


def column_residue(
    parsed: dict[str, Any], channel: int
) -> tuple[int, int]:
    moduli = [item["modulus"] for item in parsed["moduli"]]
    first_parts: list[int] = []
    second_parts: list[int] = []
    for item in parsed["moduli"]:
        residuals = [pair[channel] for pair in item["residuals"]]
        value = solve_modular_system(
            parsed["coefficients"],
            residuals,
            item["modulus"],
        )
        first_parts.append(value[0])
        second_parts.append(value[1])
    return crt(first_parts, moduli), crt(second_parts, moduli)


def candidate_vectors(
    parsed: dict[str, Any], channel: int
) -> list[tuple[int, int]]:
    modulus_product = math.prod(item["modulus"] for item in parsed["moduli"])
    residue_first, residue_second = column_residue(parsed, channel)
    bound = LIMITS["max_basis_entry"]
    candidates: list[tuple[int, int]] = []
    for first in range(-bound, bound + 1):
        if (first - residue_first) % modulus_product != 0:
            continue
        for second in range(-bound, bound + 1):
            if (second - residue_second) % modulus_product != 0:
                continue
            vector = (first, second)
            valid = True
            for item in parsed["moduli"]:
                modulus = item["modulus"]
                for coefficient, pair in zip(
                    parsed["coefficients"], item["residuals"]
                ):
                    value = (
                        coefficient[0] * vector[0]
                        + coefficient[1] * vector[1]
                    ) % modulus
                    if value != pair[channel]:
                        valid = False
                        break
                if not valid:
                    break
            if valid:
                candidates.append(vector)
    if not candidates:
        raise ContractError("no_candidate")
    return candidates


def dot(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[0] + left[1] * right[1]


def squared_norm(vector: tuple[int, int]) -> int:
    return dot(vector, vector)


def determinant(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def primitive(basis: tuple[tuple[int, int], tuple[int, int]]) -> bool:
    det = abs(determinant(*basis))
    divisor = det
    for vector in basis:
        divisor = math.gcd(divisor, abs(vector[0]))
        divisor = math.gcd(divisor, abs(vector[1]))
    return det > 0 and divisor == 1


def gauss_reduce(
    first: tuple[int, int], second: tuple[int, int]
) -> tuple[tuple[int, int], tuple[int, int]]:
    if determinant(first, second) == 0:
        raise ContractError("degenerate")
    left, right = first, second
    if (squared_norm(right), right) < (squared_norm(left), left):
        left, right = right, left
    for _ in range(32):
        denominator = squared_norm(left)
        multiplier = int(round(dot(left, right) / denominator))
        if multiplier == 0:
            break
        right = (
            right[0] - multiplier * left[0],
            right[1] - multiplier * left[1],
        )
        if (squared_norm(right), right) < (squared_norm(left), left):
            left, right = right, left
    return left, right


def canonical_pair(
    first: tuple[int, int], second: tuple[int, int]
) -> tuple[tuple[int, int], tuple[int, int]]:
    reduced_first, reduced_second = gauss_reduce(first, second)
    candidates: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for first_sign in (1, -1):
        for second_sign in (1, -1):
            oriented_first = (
                first_sign * reduced_first[0],
                first_sign * reduced_first[1],
            )
            oriented_second = (
                second_sign * reduced_second[0],
                second_sign * reduced_second[1],
            )
            if determinant(oriented_first, oriented_second) <= 0:
                continue
            if oriented_first[0] <= 0 or oriented_second[1] <= 0:
                continue
            candidates.append((oriented_first, oriented_second))
    if not candidates:
        raise ContractError("direction")
    candidates.sort(
        key=lambda basis: (
            min(squared_norm(basis[0]), squared_norm(basis[1])),
            basis[0],
            squared_norm(basis[0]),
            basis[1],
            squared_norm(basis[1]),
        )
    )
    return candidates[0]


def satisfies_all(
    basis: tuple[tuple[int, int], tuple[int, int]],
    parsed: dict[str, Any],
) -> bool:
    for item in parsed["moduli"]:
        modulus = item["modulus"]
        for coefficient, pair in zip(parsed["coefficients"], item["residuals"]):
            first = (
                coefficient[0] * basis[0][0]
                + coefficient[1] * basis[0][1]
            ) % modulus
            second = (
                coefficient[0] * basis[1][0]
                + coefficient[1] * basis[1][1]
            ) % modulus
            if first != pair[0] or second != pair[1]:
                return False
    return True


def solve(parsed: dict[str, Any]) -> dict[str, Any]:
    first_candidates = candidate_vectors(parsed, 0)
    second_candidates = candidate_vectors(parsed, 1)
    valid: list[tuple[tuple[int, int], tuple[int, int]]] = []
    for first in first_candidates:
        for second in second_candidates:
            try:
                basis = canonical_pair(first, second)
            except ContractError:
                continue
            if not primitive(basis):
                continue
            if not satisfies_all(basis, parsed):
                continue
            valid.append(basis)
    if not valid:
        raise ContractError("no_canonical_basis")

    def key(basis: tuple[tuple[int, int], tuple[int, int]]) -> tuple[Any, ...]:
        shortest = min(squared_norm(basis[0]), squared_norm(basis[1]))
        return (
            shortest,
            basis[0],
            squared_norm(basis[0]),
            basis[1],
            squared_norm(basis[1]),
        )

    basis = min(valid, key=key)
    det = determinant(*basis)
    if det <= 0 or det > LIMITS["max_determinant"]:
        raise ContractError("basis_determinant")
    fulfilled = 2 * len(parsed["coefficients"]) * len(parsed["moduli"])
    return {
        "schema_version": 1,
        "case_id": parsed["case_id"],
        "basis": [list(basis[0]), list(basis[1])],
        "check": {
            "observation_count": len(parsed["coefficients"]),
            "modulus_count": len(parsed["moduli"]),
            "satisfied_equations": fulfilled,
        },
    }


def canonical_bytes(value: dict[str, Any]) -> bytes:
    text_value = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return (text_value + "\n").encode("utf-8")


def atomic_write(path: Path, data: bytes, mode: int = 0o644) -> None:
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
    parsed = validate(load_payload(Path(args.input)))
    atomic_write(Path(args.output), canonical_bytes(solve(parsed)))


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

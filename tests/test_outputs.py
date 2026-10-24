from __future__ import annotations

import atexit
import copy
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


ENTRY = os.environ.get("YEL007_ENTRY", "/app/lattice_audit.py")
BASIS = ((1, -1), (1, 2))
COEFFICIENTS = [
    (1, 0),
    (0, 1),
    (1, 1),
    (1, -1),
    (-1, 1),
    (2, 1),
]
MODULI = [2, 3]
_WORKSPACES: list[Path] = []


def _opaque_workspace() -> Path:
    path = Path(tempfile.mkdtemp(prefix="mx-", dir="/tmp"))
    _WORKSPACES.append(path)
    return path


def _cleanup_workspaces() -> None:
    for path in _WORKSPACES:
        shutil.rmtree(path, ignore_errors=True)


atexit.register(_cleanup_workspaces)


TEST_CASE_IDS = {
    "test_recovers_canonical_basis_from_six_congruences": (
        "MODULAR-01",
        "INVERSE-26",
        "ENTRY-27",
    ),
    "test_rejects_nonminimal_positive_lift": "MODULAR-02",
    "test_near_degenerate_equal_norms_use_lex_order": (
        "MODULAR-03",
        "B1-POS-35",
        "B2-REDUCED-36",
    ),
    "test_reports_positive_primitive_determinant": (
        "MODULAR-04",
        "BASIS-ENTRY-32",
        "BASIS-DET-33",
        "DET-POS-34",
    ),
    "test_reports_two_moduli_and_equation_count": "MODULAR-05",
    "test_result_schema_is_exact": ("SCHEMA-06", "BASIS-OUTPUT-37"),
    "test_result_is_compact_sorted_json_with_one_newline": "SCHEMA-07",
    "test_success_path_keeps_console_streams_empty": "OUTPUT-08",
    "test_case_identifier_is_preserved": "OUTPUT-09",
    "test_rejects_one_modulus": "BOUNDARY-10",
    "test_rejects_three_moduli": "BOUNDARY-11",
    "test_rejects_five_observations": "BOUNDARY-12",
    "test_rejects_seven_observations": "BOUNDARY-13",
    "test_rejects_coefficient_outside_bound": "SCHEMA-14",
    "test_rejects_zero_coefficient_pair": "SCHEMA-15",
    "test_rejects_composite_modulus": "SCHEMA-16",
    "test_rejects_duplicate_modulus": "SCHEMA-17",
    "test_rejects_residual_length_mismatch": "SCHEMA-18",
    "test_rejects_inconsistent_residual_system": "MODULAR-19",
    "test_rejects_residual_outside_range": "SCHEMA-20",
    "test_rejects_wrong_schema_version": "SCHEMA-VERSION-28",
    "test_rejects_extra_modulus_field": "MODULUS-KEYS-29",
    "test_rejects_top_level_extra_key": "SCHEMA-21",
    "test_rejects_missing_limits_key": "SCHEMA-22",
    "test_invalid_input_preserves_previous_output": "RECOVERY-23",
    "test_replace_failure_preserves_directory": "RECOVERY-24",
    "test_business_result_gate": ("GATE-25", "LIMIT-ALLOW-30"),
}


def _dot(coefficient: tuple[int, int], vector: tuple[int, int]) -> int:
    return coefficient[0] * vector[0] + coefficient[1] * vector[1]


def _payload(
    case_id: str = "congruence-0712",
    basis: tuple[tuple[int, int], tuple[int, int]] = BASIS,
    moduli: list[int] | None = None,
) -> dict:
    selected_moduli = MODULI if moduli is None else moduli
    return {
        "schema_version": 1,
        "case_id": case_id,
        "coefficients": [list(value) for value in COEFFICIENTS],
        "moduli": [
            {
                "modulus": modulus,
                "residuals": [
                    [_dot(coefficient, basis[0]) % modulus, _dot(coefficient, basis[1]) % modulus]
                    for coefficient in COEFFICIENTS
                ],
            }
            for modulus in selected_moduli
        ],
        "limits": {
            "min_observations": 6,
            "max_observations": 6,
            "max_dimension": 2,
            "min_moduli": 2,
            "max_moduli": 2,
            "max_coefficient": 2,
            "max_basis_entry": 8,
            "max_determinant": 32,
        },
    }


def _invoke(
    tmp_path: Path,
    payload: dict,
    entry: str = ENTRY,
    existing_output: bytes | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    workspace = _opaque_workspace()
    input_path = workspace / "case.json"
    output_path = workspace / "basis.json"
    input_path.write_text(json.dumps(payload), encoding="utf-8")
    if existing_output is not None:
        output_path.write_bytes(existing_output)
    completed = subprocess.run(
        ["python3", entry, "--input", str(input_path), "--output", str(output_path)],
        capture_output=True,
        text=True,
        timeout=10,
        env=os.environ.copy(),
        cwd=workspace,
    )
    return completed, output_path


def _raw_invoke(tmp_path: Path, raw: bytes) -> tuple[subprocess.CompletedProcess[str], Path]:
    workspace = _opaque_workspace()
    input_path = workspace / "raw.json"
    output_path = workspace / "raw-basis.json"
    input_path.write_bytes(raw)
    completed = subprocess.run(
        ["python3", ENTRY, "--input", str(input_path), "--output", str(output_path)],
        capture_output=True,
        text=True,
        timeout=10,
        env=os.environ.copy(),
        cwd=workspace,
    )
    return completed, output_path


def _result(
    tmp_path: Path, payload: dict | None = None
) -> tuple[dict, bytes, subprocess.CompletedProcess[str]]:
    case = _payload() if payload is None else payload
    completed, output_path = _invoke(tmp_path, case)
    if not output_path.exists():
        return {}, b"", completed
    raw = output_path.read_bytes()
    return json.loads(raw), raw, completed


def test_recovers_canonical_basis_from_six_congruences(tmp_path: Path):
    """题面：「票据最外层只保留以下五个键：」"""
    # case_id: MODULAR-01, INVERSE-26, ENTRY-27
    result, _, completed = _result(tmp_path)
    assert completed.returncode == 0
    assert result["basis"] == [list(BASIS[0]), list(BASIS[1])]


def test_rejects_nonminimal_positive_lift(tmp_path: Path):
    """题面：「若有多个候选，先取最短基向量的平方范数最小值；」"""
    # case_id: MODULAR-02
    result, _, completed = _result(tmp_path)
    assert completed.returncode == 0
    assert result["basis"] != [[1, 5], [1, 2]]
    assert result["basis"][0] == [1, -1]


def test_near_degenerate_equal_norms_use_lex_order(tmp_path: Path):
    """题面：「平方范数相同时，再按 `b1`、`b2` 的字典序决胜。」"""
    # case_id: MODULAR-03, B1-POS-35, B2-REDUCED-36
    near_degenerate = ((1, -2), (2, 1))
    result, _, completed = _result(tmp_path, _payload(basis=near_degenerate))
    assert completed.returncode == 0
    assert result["basis"] == [list(near_degenerate[0]), list(near_degenerate[1])]


def test_reports_positive_primitive_determinant(tmp_path: Path):
    """题面：「恢复出的两个基向量必须使格为 primitive」"""
    # case_id: MODULAR-04, BASIS-ENTRY-32, BASIS-DET-33, DET-POS-34
    result, _, completed = _result(tmp_path)
    first = tuple(result["basis"][0])
    second = tuple(result["basis"][1])
    determinant = first[0] * second[1] - first[1] * second[0]
    divisor = determinant
    for value in (*first, *second):
        divisor = math.gcd(abs(divisor), abs(value))
    assert completed.returncode == 0
    assert determinant > 0
    assert divisor == 1


def test_reports_two_moduli_and_equation_count(tmp_path: Path):
    """题面：「`check` 只保留三个字段：」"""
    # case_id: MODULAR-05
    result, _, completed = _result(tmp_path)
    assert completed.returncode == 0
    assert result["check"] == {
        "observation_count": 6,
        "modulus_count": 2,
        "satisfied_equations": 24,
    }


def test_result_schema_is_exact(tmp_path: Path):
    """题面：「结果最外层只保留四个键：」"""
    # case_id: SCHEMA-06, BASIS-OUTPUT-37
    result, _, completed = _result(tmp_path)
    assert completed.returncode == 0
    assert set(result) == {"schema_version", "case_id", "basis", "check"}
    assert all(len(vector) == 2 for vector in result["basis"])


def test_result_is_compact_sorted_json_with_one_newline(tmp_path: Path):
    """题面：「结果使用无空格 JSON，键名升序，文件只以一个换行收尾。」"""
    # case_id: SCHEMA-07
    _, raw, completed = _result(tmp_path)
    assert completed.returncode == 0
    assert b", " not in raw
    assert raw.endswith(b"\n")
    assert not raw.endswith(b"\n\n")
    assert list(json.loads(raw)) == sorted(json.loads(raw))


def test_success_path_keeps_console_streams_empty(tmp_path: Path):
    """题面：「成功完成时 `stdout` 与 `stderr` 都必须为空字符串。」"""
    # case_id: OUTPUT-08
    _, _, completed = _result(tmp_path)
    assert completed.returncode == 0
    assert completed.stdout == ""
    assert completed.stderr == ""


def test_case_identifier_is_preserved(tmp_path: Path):
    """题面：「`schema_version` 和 `case_id` 原样沿用输入。」"""
    # case_id: OUTPUT-09
    result, _, completed = _result(tmp_path, _payload("congruence-night-7"))
    assert completed.returncode == 0
    assert result["case_id"] == "congruence-night-7"


def test_rejects_one_modulus(tmp_path: Path):
    """题面：「`moduli` 恰好有 2 个模数视角。」"""
    # case_id: BOUNDARY-10
    payload = copy.deepcopy(_payload())
    payload["moduli"].pop()
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_three_moduli(tmp_path: Path):
    """题面：「`moduli` 恰好有 2 个模数视角。」"""
    # case_id: BOUNDARY-11
    payload = copy.deepcopy(_payload())
    extra = copy.deepcopy(payload["moduli"][0])
    extra["modulus"] = 5
    extra["residuals"] = [[0, 0] for _ in extra["residuals"]]
    payload["moduli"].append(extra)
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_five_observations(tmp_path: Path):
    """题面：「`coefficients` 恰好有 6 个有序的二维整数系数」"""
    # case_id: BOUNDARY-12
    payload = copy.deepcopy(_payload())
    payload["coefficients"].pop()
    for item in payload["moduli"]:
        item["residuals"].pop()
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_seven_observations(tmp_path: Path):
    """题面：「`coefficients` 恰好有 6 个有序的二维整数系数」"""
    # case_id: BOUNDARY-13
    payload = copy.deepcopy(_payload())
    payload["coefficients"].append([1, 2])
    for item in payload["moduli"]:
        item["residuals"].append([0, 0])
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_coefficient_outside_bound(tmp_path: Path):
    """题面：「两个分量都必须位于闭区间 `[-2,2]`，且不能同时为零。」"""
    # case_id: SCHEMA-14
    payload = copy.deepcopy(_payload())
    payload["coefficients"][0] = [3, 0]
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_zero_coefficient_pair(tmp_path: Path):
    """题面：「两个分量都必须位于闭区间 `[-2,2]`，且不能同时为零。」"""
    # case_id: SCHEMA-15
    payload = copy.deepcopy(_payload())
    payload["coefficients"][0] = [0, 0]
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_composite_modulus(tmp_path: Path):
    """题面：「`modulus` 必须是互不相同的质数。」"""
    # case_id: SCHEMA-16
    payload = copy.deepcopy(_payload())
    payload["moduli"][0]["modulus"] = 4
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_duplicate_modulus(tmp_path: Path):
    """题面：「`modulus` 必须是互不相同的质数。」"""
    # case_id: SCHEMA-17
    payload = copy.deepcopy(_payload())
    payload["moduli"][1]["modulus"] = payload["moduli"][0]["modulus"]
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_residual_length_mismatch(tmp_path: Path):
    """题面：「`residuals` 的长度必须等于系数数」"""
    # case_id: SCHEMA-18
    payload = copy.deepcopy(_payload())
    payload["moduli"][0]["residuals"].pop()
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_inconsistent_residual_system(tmp_path: Path):
    """题面：「坏票据、不可解同余系统、越界基和方向规范不唯一统一返回状态码 `3`」"""
    # case_id: MODULAR-19
    payload = copy.deepcopy(_payload())
    payload["moduli"][0]["residuals"][0][0] ^= 1
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_residual_outside_range(tmp_path: Path):
    """题面：「第 `i` 项恰好是两个整数 `[left, right]`，都在 `[0, modulus-1]` 内。」"""
    # case_id: SCHEMA-20
    payload = copy.deepcopy(_payload())
    payload["moduli"][0]["residuals"][0][0] = 2
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_wrong_schema_version(tmp_path: Path):
    """题面：「`schema_version` 必须是整数 `1`。」"""
    # case_id: SCHEMA-VERSION-28
    payload = copy.deepcopy(_payload())
    payload["schema_version"] = 2
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_extra_modulus_field(tmp_path: Path):
    """题面：「每个对象的键恰好为：」"""
    # case_id: MODULUS-KEYS-29
    payload = copy.deepcopy(_payload())
    payload["moduli"][0]["quality"] = 1
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_top_level_extra_key(tmp_path: Path):
    """题面：「票据最外层只保留以下五个键：」"""
    # case_id: SCHEMA-21
    payload = copy.deepcopy(_payload())
    payload["operator"] = "audit"
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_rejects_missing_limits_key(tmp_path: Path):
    """题面：「`limits` 只保留八个规模与界值字段：」"""
    # case_id: SCHEMA-22
    payload = copy.deepcopy(_payload())
    payload.pop("limits")
    completed, output_path = _invoke(tmp_path, payload)
    assert completed.returncode == 3
    assert not output_path.exists()


def test_invalid_input_preserves_previous_output(tmp_path: Path):
    """题面：「坏票据、不可解同余系统、越界基和方向规范不唯一统一返回状态码 `3`」"""
    # case_id: RECOVERY-23
    payload = copy.deepcopy(_payload())
    payload["coefficients"][0] = [3, 0]
    completed, output_path = _invoke(
        tmp_path,
        payload,
        existing_output=b"old-basis\n",
    )
    assert completed.returncode == 3
    assert output_path.read_bytes() == b"old-basis\n"


def test_replace_failure_preserves_directory(tmp_path: Path):
    """题面：「结果经同目录临时文件安装；若最终换名失败，既有文件按旧内容保留。」"""
    # case_id: RECOVERY-24
    workspace = _opaque_workspace()
    input_path = workspace / "case.json"
    output_path = workspace / "basis.json"
    input_path.write_text(json.dumps(_payload()), encoding="utf-8")
    output_path.mkdir()
    completed = subprocess.run(
        ["python3", ENTRY, "--input", str(input_path), "--output", str(output_path)],
        capture_output=True,
        text=True,
        timeout=10,
        env=os.environ.copy(),
        cwd=workspace,
    )
    assert completed.returncode == 3
    assert output_path.is_dir()


def test_business_result_gate(tmp_path: Path):
    """题面：「审计程序要从观测恢复格，再输出唯一方向规范基。」"""
    # case_id: GATE-25, LIMIT-ALLOW-30
    result, _, completed = _result(tmp_path)
    assert completed.returncode == 0
    assert result["basis"] == [list(BASIS[0]), list(BASIS[1])]
    assert result["check"]["satisfied_equations"] == 24


def test_candidate_paths_do_not_leak_pytest_identity(tmp_path: Path):
    """候选可见路径不得包含测试名、pytest 临时目录或节点标识。"""
    completed, output_path = _invoke(tmp_path, _payload())
    assert completed.returncode == 0
    visible = output_path.as_posix()
    assert "test_" not in visible
    assert "pytest" not in visible
    assert json.loads(output_path.read_text())["case_id"] == "congruence-0712"

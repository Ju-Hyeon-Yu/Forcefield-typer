# 역장별 타이퍼(DREIDING·GAFF2·OPLS-AA) 목록과 공통 예외
"""forcefield-typer — SMILES(+mol2)에서 moltemplate 단일 분자 .lt 를 만든다.

역장 파일(dreiding.lt, gaff2.lt, oplsaa2024.lt)은 moltemplate 배포본을 그대로 읽기만 한다.
"""
__version__ = "0.1.0"

FORCE_FIELDS = ("dreiding", "gaff2", "oplsaa")


class TypingError(Exception):
    """역장 규칙으로 라벨을 정할 수 없을 때."""


def get_typer(name, ff_file=None):
    """역장 이름 → 타이퍼 객체 (assign(mol) → (라벨, 경고))."""
    if name == "dreiding":
        from .dreiding import DreidingTyper
        return DreidingTyper(ff_file)
    if name == "gaff2":
        from .gaff2 import Gaff2Typer
        return Gaff2Typer(ff_file)
    if name == "oplsaa":
        from .oplsaa import OplsTyper
        return OplsTyper(ff_file)
    raise ValueError(f"알 수 없는 역장: {name} (가능: {', '.join(FORCE_FIELDS)})")


def moltemplate_ff_path(filename):
    """설치된 moltemplate 의 force_fields/<filename> 경로."""
    import os
    import moltemplate
    return os.path.join(os.path.dirname(moltemplate.__file__), "force_fields", filename)

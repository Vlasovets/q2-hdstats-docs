from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Any


def sha256(path: str | Path) -> str | None:
    p = Path(path)
    if not p.exists() or not p.is_file():
        return None
    h = hashlib.sha256()
    with p.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def file_record(path: str | Path, project_root: str | Path | None = None) -> dict[str, Any]:
    p = Path(path)
    if project_root is not None:
        root = Path(project_root)
        try:
            label = str(p.resolve().relative_to(root.resolve()))
        except ValueError:
            label = str(p)
    else:
        label = str(p)
    return {
        'path': label,
        'absolute_path': str(p),
        'exists': p.exists(),
        'size_bytes': p.stat().st_size if p.exists() and p.is_file() else None,
        'sha256': sha256(p),
    }


def git_info(project_root: str | Path) -> dict[str, Any]:
    root = Path(project_root)
    try:
        head = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True, stderr=subprocess.DEVNULL).strip()
        status = subprocess.check_output(['git', '-C', str(root), 'status', '--short'], text=True, stderr=subprocess.DEVNULL)
    except Exception:
        return {'head': None, 'status_short': None}
    return {'head': head, 'status_short': status.splitlines()}


def write_provenance(
    output_base: str | Path,
    *,
    description: str,
    project_root: str | Path,
    script: str | Path,
    inputs: Iterable[str | Path],
    outputs: Iterable[str | Path],
    parameters: Mapping[str, Any] | None = None,
    notes: Iterable[str] | None = None,
) -> Path:
    project_root = Path(project_root)
    script = Path(script)
    payload = {
        'description': description,
        'generated_at_utc': datetime.now(timezone.utc).isoformat(),
        'project_root': str(project_root),
        'script': file_record(script, project_root),
        'command': f"python3 {script}",
        'parameters': dict(parameters or {}),
        'inputs': [file_record(p, project_root) for p in inputs],
        'outputs': [file_record(p, project_root) for p in outputs],
        'environment': {
            'python': platform.python_version(),
            'platform': platform.platform(),
            'cwd': os.getcwd(),
        },
        'git': git_info(project_root),
        'notes': list(notes or []),
    }
    out = Path(output_base).with_suffix('.provenance.json')
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    return out

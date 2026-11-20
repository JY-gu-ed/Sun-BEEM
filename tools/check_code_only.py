"""Validate that the release contains source code and documentation only."""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SPECIAL = {'.gitignore', '.gitattributes', '.env.example', 'requirements-R.txt', 'CITATION.cff'}
EXTENSIONS = {'.R', '.py', '.md', '.yml'}
RULES = {
    'credential token': re.compile(r'(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9]{24,})'),
    'embedded participant IDs': re.compile(r'\beid\s*(?:=|<-)\s*c\(\s*["\']?\d', re.I),
    'literal participant filter': re.compile(r'\beid\s*==\s*["\']?\d', re.I),
    'literal credential value': re.compile(r'(?:password|passwd|api[_-]?key|access[_-]?token|secret)\s*(?:<-|=|:)\s*["\'][^"\']+["\']', re.I),
}
failures = []
files = []
for path in sorted(ROOT.rglob('*')):
    relative = path.relative_to(ROOT)
    if '.git' in relative.parts:
        continue
    if path.is_symlink():
        failures.append(f'{relative}: symbolic links are not permitted')
        continue
    if not path.is_file():
        continue
    files.append(path)
    if path.name not in SPECIAL and path.suffix not in EXTENSIONS:
        failures.append(f'{relative}: prohibited file type')
        continue
    raw = path.read_bytes()
    if b'\x00' in raw:
        failures.append(f'{relative}: binary content')
        continue
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        failures.append(f'{relative}: not UTF-8 text')
        continue
    if len(raw) > 300_000:
        failures.append(f'{relative}: unexpectedly large text file; review before release')
    for number, line in enumerate(text.splitlines(), 1):
        for label, rule in RULES.items():
            if rule.search(line):
                failures.append(f'{relative}:{number}: {label}')
    if path.suffix == '.R' and re.search(r'/home/[^\s"\']+|[A-Z]:[/\\]Users[/\\]', text, re.I):
        failures.append(f'{relative}: private workstation path')

if failures:
    print('Release check failed. Contents are not printed:')
    print('\n'.join(failures))
    sys.exit(1)
print(f'PASS: {len(files)} UTF-8 source/documentation files; no prohibited data files, embedded participant IDs or credential patterns detected.')

"""Check the supplement/web contract without recomputing statistical models."""
from pathlib import Path
import hashlib
import json
import re
import shlex
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qiime_commands(text):
    result = []
    for block in re.findall(r'^```bash\n(.*?)\n```', text, re.M | re.S):
        for line in block.replace('\\\n', ' ').splitlines():
            if line.strip().startswith('qiime '):
                result.append(shlex.split(line, comments=True))
    return result


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else '.').resolve()
    bundle = root / 'analysis/supplement-v2'
    manifest = json.loads((bundle / 'conversion-manifest.json').read_text())
    source = bundle / 'source/supplementary.tex'
    assert sha(source) == manifest['source_sha256'], 'Supplement changed: regenerate the web pages.'
    blocks = []
    pages = []
    editorial_path = bundle / 'editorial-manifest.json'
    editorial = json.loads(editorial_path.read_text()) if editorial_path.exists() else {}
    for page in manifest['sections']:
        path = root / 'docs' / (page['page'] + '.md')
        text = path.read_text()
        pages.append(text)
        if page['page'] in editorial:
            revision = editorial[page['page']]
            baseline = bundle / revision['baseline']
            assert sha(baseline) == page['sha256'], f'Supplement baseline drift: {path}'
            assert sha(path) == revision['published_sha256'], f'Unreviewed tutorial edit: {path}'
            published_commands = qiime_commands(text)
            replacements = revision.get('command_filename_replacements', {})
            for command in qiime_commands(baseline.read_text()):
                expected = [replacements.get(token, token) for token in command]
                assert expected in published_commands, f'Source analysis command missing or changed: {expected}'
            blocks += re.findall(r'^```bash\n(.*?)\n```', baseline.read_text(), flags=re.M | re.S)
        else:
            assert sha(path) == page['sha256'], f'Generated page drift: {path}'
            blocks += re.findall(r'^```bash\n(.*?)\n```', text, flags=re.M | re.S)
    source_blocks = re.findall(r'\\begin\{(?:PseudoCode|CodeInput|CodeOutput)\}\n(.*?)\n\\end\{(?:PseudoCode|CodeInput|CodeOutput)\}', source.read_text(), flags=re.S)
    assert blocks == [b.strip('\n') for b in source_blocks], 'A command block differs from the supplement.'
    assert len(blocks) == 30 and len(manifest['figures']) == 26 and len(manifest['tables']) == 10
    for figure in manifest['figures']:
        pdf = root / 'docs/images/tutorial-v2' / figure['pdf']
        assert sha(pdf) == figure['sha256'] == sha(bundle / 'source/figs' / pdf.name)
        assert sha(pdf.with_suffix('.png')) == figure['png_sha256']
        assert pdf.with_suffix('.provenance.json').is_file()
    for path, expected in json.loads((bundle / 'branding-manifest.json').read_text()).items():
        assert sha(root / path) == expected, 'Logo changed: ' + path
    source_keys = set()
    for item in re.findall(r'\\cite(?:p|t|yearpar)?(?:\[[^]]*\])*\{([^}]+)\}', source.read_text()):
        source_keys.update(k.strip() for k in item.split(','))
    rendered = '\n'.join(pages)
    rendered_keys = set()
    for item in re.findall(r'\{cite(?::\w+)?\}`([^`]+)`', rendered):
        rendered_keys.update(k.strip() for k in item.split(','))
    rendered_keys.update(re.findall(r'href="../references.html#([^"]+)"', rendered))
    assert source_keys == rendered_keys, f'Citation mismatch: {source_keys ^ rendered_keys}'
    uncommented = re.sub(r'(?<!\\)%[^\n]*', '', source.read_text())
    source_labels = set(re.findall(r'\\label\{([^}]+)\}', uncommented))
    rendered_labels = set(re.findall(r'^\(([^\n]+)\)=$|^:(?:label|name): ([^\n]+)$', rendered, re.M))
    targets = {a or b for a, b in rendered_labels}
    assert not (source_labels - targets), f'Lost labels: {source_labels - targets}'
    policy = bundle / 'command-harmonization'
    if policy.exists():
        from minimize_qiime_commands import simplify_text, resolved
        defaults_path = policy / 'qiime-interface-defaults.json'
        defaults = json.loads(defaults_path.read_text())
        report = json.loads((policy / 'report.json').read_text())
        assert sha(defaults_path) == report['default_snapshot_sha256'], 'Default snapshot changed'
        texts = [(source.read_text(), True)] + [(page, False) for page in pages]
        for text, latex in texts:
            _, commands = simplify_text(text, defaults, latex)
            assert not any(cmd['removed'] for cmd in commands), 'Redundant default parameter in displayed command'
        for file in report['files']:
            for command in file['commands']:
                assert resolved(command['before'], defaults) == resolved(command['after'], defaults), 'Effective command settings changed'
    if 'tutorial/s2-01' in editorial:
        import csv
        import io
        import zipfile
        with zipfile.ZipFile(root / 'docs/downloads/atacama-top15-inputs.zip') as archive:
            labels = csv.DictReader(io.StringIO(archive.read('top15-asv-label-map.tsv').decode()), delimiter='\t')
            lowdim = (root / 'docs/tutorial/s2-01.md').read_text()
            for row in labels:
                assert row['feature_id'] in lowdim, 'Missing ASV feature identifier: ' + row['asv_label']
            inputs = json.loads(archive.read('input-manifest.json'))
            for item in inputs['files']:
                assert hashlib.sha256(archive.read(item['file'])).hexdigest() == item['sha256']
        with zipfile.ZipFile(root / 'docs/downloads/atacama-top15-plotting.zip') as archive:
            for member in archive.namelist():
                relative = Path(member).relative_to('plotting')
                bundled = bundle / 'portable-plotting/top15' / relative
                assert archive.read(member) == bundled.read_bytes(), 'Download/script mismatch: ' + member
    print(f'PASS: {len(pages)} pages; source commands preserved; {len(editorial)} reviewed tutorial adaptations; 26 unchanged PDFs; 10 tables; {len(source_keys)} citation keys; {len(source_labels)} source labels; both original logos; minimal-command policy checked.')


if __name__ == '__main__':
    main()

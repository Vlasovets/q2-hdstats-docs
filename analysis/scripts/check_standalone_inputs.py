"""Check public input ZIP integrity without a QIIME installation."""
from pathlib import Path
import hashlib
import json
import sys
import zipfile


def main():
    repo = Path(sys.argv[1] if len(sys.argv) > 1 else '.').resolve()
    records = repo / 'analysis/supplement-v2/standalone-inputs'
    bundles = json.loads((records / 'bundle-manifest.json').read_text())
    count = 0
    for bundle, spec in bundles.items():
        path = repo / 'docs/downloads' / (bundle + '.zip')
        assert hashlib.sha256(path.read_bytes()).hexdigest() == spec['zip_sha256'], bundle
        with zipfile.ZipFile(path) as archive:
            manifest = json.loads(archive.read('input-manifest.json'))
            assert manifest == json.loads((records / (bundle + '-manifest.json')).read_text())
            assert manifest['dimensions'] == spec['dimensions']
            assert len(manifest['files']) == spec['files']
            for item in manifest['files']:
                name = item['file']
                assert not Path(name).is_absolute() and '..' not in Path(name).parts
                content = archive.read(name)
                assert hashlib.sha256(content).hexdigest() == item['sha256'], name
                assert len(content) == item['bytes'], name
                if name.endswith('.py'):
                    compile(content, name, 'exec')
                if name.endswith('.qza'):
                    import io
                    with zipfile.ZipFile(io.BytesIO(content)) as artifact:
                        metadata = next(n for n in artifact.namelist() if n.count('/') == 1 and n.endswith('/metadata.yaml'))
                        artifact_type = artifact.read(metadata).decode()
                        assert 'FeatureTable[Frequency]' in artifact_type or 'FeatureData[Taxonomy]' in artifact_type, name
                count += 1
    validation = json.loads((records / 'input-validation.json').read_text())
    assert validation and all(item['passed'] for item in validation)
    print(f'PASS: {len(bundles)} public input bundles; {count} hashed files; frequency/taxonomy artifacts only; recorded preparation checks pass.')


if __name__ == '__main__':
    main()

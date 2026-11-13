"""Validate local links and image targets in the supplement-derived HTML."""
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    def __init__(self, path):
        super().__init__(convert_charrefs=True)
        self.ids = set()
        self.links = []
        self.feed(path.read_text())

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if 'id' in attrs:
            self.ids.add(attrs['id'])
        if tag == 'a' and attrs.get('href'):
            self.links.append(attrs['href'])
        if tag == 'img' and attrs.get('src'):
            self.links.append(attrs['src'])


def main():
    repo = Path(sys.argv[1] if len(sys.argv) > 1 else '.').resolve()
    root = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else repo / 'docs/_build/html'
    manifest = json.loads((repo / 'analysis/supplement-v2/conversion-manifest.json').read_text())
    pages = [root / (p['page'] + '.html') for p in manifest['sections']]
    pages += [root / (name + '.html') for name in ('intro', 'getting-started', 'references', 'reproducibility')]
    cache = {}
    errors = []
    count = 0
    for path in pages:
        cache[path] = Page(path)
        for target in cache[path].links:
            url = urlsplit(target)
            if url.scheme or url.netloc or not (url.path or url.fragment):
                continue
            destination = (path.parent / unquote(url.path)).resolve() if url.path else path
            count += 1
            if not destination.exists():
                errors.append((str(path), target, 'missing file'))
                continue
            if destination.suffix == '.html' and url.fragment:
                if destination not in cache:
                    cache[destination] = Page(destination)
                if unquote(url.fragment) not in cache[destination].ids:
                    errors.append((str(path), target, 'missing anchor'))
    if errors:
        raise SystemExit(json.dumps(errors, indent=2))
    print(f'PASS: {len(pages)} built pages; {count} local links/images; no missing files or anchors.')


if __name__ == '__main__':
    main()

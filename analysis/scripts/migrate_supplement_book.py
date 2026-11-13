"""Convert the authoritative LaTeX supplement to MyST through Pandoc's AST.

No models are fitted. Code blocks and source figures are preserved. Run against
a staged repository first, then review its conversion report and strict build.
"""
import argparse
import collections
import copy
import hashlib
import html
import json
from pathlib import Path
import re
import shutil
import subprocess


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def balanced(text, start):
    assert text[start] == '{'
    depth = 1
    pos = start + 1
    while depth:
        if text[pos] == '{' and text[pos - 1] != '\\':
            depth += 1
        elif text[pos] == '}' and text[pos - 1] != '\\':
            depth -= 1
        pos += 1
    return text[start + 1:pos - 1], pos


def nodes(value):
    if isinstance(value, dict):
        if 't' in value:
            yield value
        for child in value.values():
            yield from nodes(child)
    elif isinstance(value, list):
        for child in value:
            yield from nodes(child)


def raw(text, block=False):
    return {'t': 'RawBlock' if block else 'RawInline', 'c': ['markdown', text]}


class Converter:
    def __init__(self, root, source):
        self.root, self.source = root, source
        self.docs = root / 'docs'
        self.dest = self.docs / 'tutorial'
        self.assets = self.docs / 'images/tutorial-v2'
        self.record = {'source_sha256': sha(source), 'sections': [], 'figures': [],
                       'commands': [], 'labels': [], 'tables': []}
        self.record['pandoc_version'] = subprocess.run(['pandoc', '--version'], text=True, capture_output=True, check=True).stdout.splitlines()[0]
        self.record['poppler_version'] = subprocess.run(['pdftoppm', '-v'], text=True, capture_output=True, check=True).stderr.splitlines()[0]
        self.code = []
        self.table_captions = []
        self.figure_count = 0
        self.table_count = 0
        self.known_labels = set()
        self.page_for_label = {}
        self.number_for_label = {}

    def preprocess(self, text):
        # These are listing environments, not LaTeX prose. Protect them before
        # parsing so shell backslashes and spaces remain byte-for-byte intact.
        def code(match):
            body = match.group(2).strip('\n')
            self.code.append(body)
            return '\n\\begin{verbatim}\n' + body + '\n\\end{verbatim}\n'
        text = re.sub(r'\\begin\{(PseudoCode|CodeInput|CodeOutput)\}(.*?)\\end\{\1\}',
                      code, text, flags=re.S)
        text = text.replace(r'\newcommand{\code}[1]{\allowbreak\texttt{\hyphenchar\font=`\-\relax #1}}',
                            r'\newcommand{\code}[1]{\texttt{#1}}')
        text = text.replace(r'\newcommand{\asvid}[1]{\raisebox{0.15ex}{\code{#1}}}',
                            r'\newcommand{\asvid}[1]{\texttt{#1}}')
        # Pandoc supports table contents/spans but not every physical TeX
        # column specification. Normalize layout only, preserving all cells.
        def table(match):
            block = match.group(0)
            spec = re.search(r'\\begin\{(tabularx|tabular|longtable)\}', block)
            if not spec:
                return block
            pos = spec.end()
            if spec.group(1) == 'tabularx':
                _, pos = balanced(block, pos)
            _, end = balanced(block, pos)
            tail = block[end:]
            header = re.search(r'\\textbf\{[^\n]+', tail)
            if not header:
                raise ValueError('Table has no recognizable header')
            columns = header.group(0).split(r'\\')[0].count('&') + 1
            begin = '\\begin{' + ('longtable' if spec.group(1) == 'longtable' else 'tabular') + '}{' + 'l' * columns + '}'
            block = block[:spec.start()] + begin + block[end:]
            block = block.replace(r'\end{tabularx}', r'\end{tabular}')
            block = block.replace('@{}', '')
            # Preserve Pandoc's recognized longtable header/footer markers.
            return block
        text = re.sub(r'\\begin\{table\}.*?\\end\{table\}|\\begin\{longtable\}.*?\\end\{longtable\}',
                      table, text, flags=re.S)
        return text

    def pandoc(self, blocks, target='markdown'):
        doc = {'pandoc-api-version': self.api, 'meta': {}, 'blocks': blocks}
        command = ['pandoc', '-f', 'json', '-t', target, '--wrap=none']
        if target == 'markdown':
            command += ['--markdown-headings=atx']
        return subprocess.run(command, input=json.dumps(doc), capture_output=True,
                              text=True, check=True).stdout.strip()

    def inline(self, values):
        return self.pandoc([{'t': 'Plain', 'c': self.transform(values)}])

    def transform(self, obj):
        if isinstance(obj, list):
            return [self.transform(x) for x in obj]
        if not isinstance(obj, dict):
            return obj
        t, c = obj.get('t'), obj.get('c')
        if t == 'Cite':
            keys = ','.join(x['citationId'] for x in c[0])
            role = 'cite:t' if any(x['citationMode']['t'] == 'AuthorInText' for x in c[0]) else 'cite:p'
            if any(x['citationMode']['t'] == 'SuppressAuthor' for x in c[0]):
                role = 'cite:yearpar'
            return raw('{' + role + '}`' + keys + '`')
        if t == 'Span' and c[0][0]:
            if c[1]:
                raise ValueError('Unexpected nonempty labeled span')
            return raw('')
        if t == 'Link' and c[2][0].startswith('#'):
            label = c[2][0][1:]
            if label not in self.known_labels:
                raise ValueError('Unresolved LaTeX reference: ' + label)
            role = 'eq' if label.startswith('eq:') else 'ref'
            if role == 'ref':
                kind, number = self.number_for_label.get(label, ('Section', ''))
                reference_type = dict(c[0][2]).get('reference-type', '')
                title = (kind + ' ' if reference_type == 'autoref' else '') + number
                return raw('{ref}`' + title + ' <' + label + '>`')
            return raw('{' + role + '}`' + label + '`')
        if t == 'Math' and c[0]['t'] == 'DisplayMath':
            labels = re.findall(r'\\label\{([^}]+)\}', c[1])
            if len(labels) > 1:
                inner = c[1].removeprefix(r'\begin{aligned}').removesuffix(r'\end{aligned}').strip()
                rows = re.split(r'\\\\(?:\[[^]]*\])?', inner)
                if len(rows) != len(labels) or any(len(re.findall(r'\\label\{([^}]+)\}', row)) != 1 for row in rows):
                    raise ValueError('Unsupported multiply labeled equation block')
                return raw('\n\n'.join(self.transform({'t': 'Math', 'c': [c[0], r'\begin{aligned}' + row + r'\end{aligned}']})['c'][1] for row in rows))
            formula = re.sub(r'\\label\{[^}]+\}', '', c[1]).strip()
            directive = '\n\n```{math}\n'
            if labels:
                directive += ':label: ' + labels[0] + '\n'
            directive += '\n' + formula + '\n```\n\n'
            return raw(directive)
        if t == 'RawInline':
            if c[1] in (r'\allowbreak', r'\relax'):
                return raw('')
            raise ValueError('Unconverted inline LaTeX: ' + str(c))
        return {key: self.transform(value) for key, value in obj.items()}

    def table(self, node, label=''):
        table = copy.deepcopy(node)
        if table['t'] == 'Div':
            label = table['c'][0][0] or label
            children = table['c'][1]
            if len(children) != 1 or children[0]['t'] != 'Table':
                return self.blocks(children)
            table = children[0]
        c = table['c']
        label = label or c[0][0]
        self.table_count += 1
        caption = self.inline(c[1][1][0]['c']) if c[1][1] else ''
        # HTML preserves multicolumn spans; Pandoc's MathJax spans retain all
        # mathematics. Table citations link to stable bibliography key anchors.
        def html_nodes(obj):
            if isinstance(obj, list):
                return [html_nodes(x) for x in obj]
            if not isinstance(obj, dict):
                return obj
            if obj.get('t') == 'Cite':
                keys = [v['citationId'] for v in obj['c'][0]]
                return {'t': 'RawInline', 'c': ['html', ', '.join(
                    '<a href="../references.html#' + html.escape(k) + '">' + html.escape(k) + '</a>' for k in keys)]}
            if obj.get('t') == 'Link' and obj['c'][2][0].startswith('#'):
                label = obj['c'][2][0][1:]
                page = self.page_for_label.get(label)
                if not page:
                    raise ValueError('Unresolved table reference: ' + label)
                obj = copy.deepcopy(obj)
                anchor = ('equation-' if label.startswith('eq:') else '') + label.replace(':', '-').lower()
                obj['c'][2][0] = page + '.html#' + anchor
                obj['c'][1] = [{'t': 'Str', 'c': label}]
            return {k: html_nodes(v) for k, v in obj.items()}
        c[1] = [None, []]
        table = html_nodes(table)
        doc = {'pandoc-api-version': self.api, 'meta': {}, 'blocks': [table]}
        result = subprocess.run(['pandoc', '-f', 'json', '-t', 'html', '--mathjax'],
                                input=json.dumps(doc), text=True, capture_output=True, check=True).stdout
        self.record['tables'].append({'number': self.table_count, 'label': label})
        anchor = '(' + label + ')=\n' if label else ''
        return anchor + '**Table S' + str(self.table_count) + '.** ' + caption + '\n\n' + result

    def blocks(self, blocks):
        output = []
        for block in blocks:
            t, c = block['t'], block.get('c')
            if t in ('Para', 'Plain'):
                for node in nodes(block):
                    if node['t'] == 'Span' and node['c'][0][0]:
                        output.append('(' + node['c'][0][0] + ')=')
            if t == 'Div':
                if any(x['t'] == 'Table' for x in c[1]):
                    output.append(self.table(block))
                else:
                    if c[0][0]:
                        output.append('(' + c[0][0] + ')=')
                    output.append(self.blocks(c[1]))
            elif t == 'Figure':
                images = [x for x in nodes(block) if x['t'] == 'Image']
                if len(images) != 1:
                    raise ValueError('Expected one source PDF per figure')
                image = images[0]['c'][2][0]
                pdf = self.source.parent / image
                stem = pdf.stem
                shutil.copy2(pdf, self.assets / pdf.name)
                provenance = pdf.with_suffix('.provenance.json')
                if not provenance.exists():
                    raise ValueError('Missing audited figure provenance: ' + str(pdf))
                shutil.copy2(provenance, self.assets / provenance.name)
                png = self.assets / (stem + '.png')
                subprocess.run(['pdftoppm', '-singlefile', '-scale-to', '2400', '-png', str(pdf), str(png.with_suffix(''))], check=True)
                self.figure_count += 1
                caption = self.blocks(c[1][1])
                label = c[0][0] or images[0]['c'][0][0]
                output.append('```{figure} /images/tutorial-v2/' + png.name + '\n:name: ' + label + '\n:width: 100%\n\nFigure S' + str(self.figure_count) + '. ' + caption + '\n```\n\n[PDF](/images/tutorial-v2/' + pdf.name + ') | [Provenance](/images/tutorial-v2/' + provenance.name + ')')
                self.record['figures'].append({'number': self.figure_count, 'label': label, 'pdf': pdf.name, 'sha256': sha(pdf), 'png_sha256': sha(png)})
            elif t == 'Table':
                output.append(self.table(block))
            elif t == 'Header':
                level, attrs, title = c
                label = attrs[0]
                if level == 4:
                    output.append('**' + self.inline(title) + '**')
                else:
                    output.append(('(' + label + ')=\n' if label else '') + '#' * max(1, level - self.page_level + 1) + ' ' + self.inline(title))
            elif t == 'CodeBlock':
                code = c[1]
                self.record['commands'].append({'sha256': hashlib.sha256(code.encode()).hexdigest(), 'text': code})
                output.append('```bash\n' + code + '\n```')
            elif t == 'RawBlock':
                if c[1].strip() not in (r'\FloatBarrier', r'\clearpage', r'\newpage'):
                    raise ValueError('Unconverted block LaTeX: ' + str(c))
            else:
                output.append(self.pandoc([self.transform(block)]))
        return '\n\n'.join(x for x in output if x.strip())

    def run(self):
        editorial = self.root / 'analysis/supplement-v2/editorial-manifest.json'
        if editorial.exists():
            raise RuntimeError('This book has reviewed tutorial adaptations. Convert into a clean staging copy and review the changes; do not overwrite the edited pages.')
        self.dest.mkdir(parents=True, exist_ok=True)
        self.assets.mkdir(parents=True, exist_ok=True)
        shutil.copy2(self.source.with_suffix('.bib'), self.docs / 'supplementary.bib')
        source = self.source.read_text()
        preprocessed = self.preprocess(source)
        result = subprocess.run(['pandoc', '-f', 'latex', '-t', 'json'], input=preprocessed,
                                text=True, capture_output=True, check=True)
        doc = json.loads(result.stdout)
        self.api = doc['pandoc-api-version']
        blocks = doc['blocks']
        first = next(i for i, b in enumerate(blocks) if b['t'] == 'Header' and b['c'][0] == 1)
        blocks = blocks[first:]
        pages = []
        section = 0
        subsection = 0
        for block in blocks:
            if block['t'] == 'Header' and block['c'][0] <= 2:
                level = block['c'][0]
                if level == 1:
                    section += 1
                    subsection = 0
                else:
                    subsection += 1
                name = f's{section}' + (f'-{subsection:02d}' if subsection else '')
                pages.append({'name': name, 'level': level, 'blocks': [], 'section': section})
            pages[-1]['blocks'].append(block)
        for page in pages:
            for node in nodes(page['blocks']):
                attrs = node['c'][1] if node['t'] == 'Header' else node.get('c', [None])[0]
                if node['t'] in ('Header', 'Figure', 'Div', 'Span') and attrs[0]:
                    self.known_labels.add(attrs[0])
                    self.page_for_label[attrs[0]] = page['name']
                if node['t'] == 'Math':
                    for label in re.findall(r'\\label\{([^}]+)\}', node['c'][1]):
                        self.known_labels.add(label)
                        self.page_for_label[label] = page['name']
        figure_number = table_number = 0
        for page in pages:
            for node in nodes(page['blocks']):
                if node['t'] == 'Figure':
                    figure_number += 1
                    self.number_for_label[node['c'][0][0]] = ('Figure', f'S{figure_number}')
                elif node['t'] == 'Div' and any(x['t'] == 'Table' for x in node['c'][1]):
                    table_number += 1
                    self.number_for_label[node['c'][0][0]] = ('Table', f'S{table_number}')
                elif node['t'] == 'Header':
                    label = node['c'][1][0]
                    number = 'S' + str(page['section'])
                    if page['level'] == 2:
                        number += '.' + str(int(page['name'].split('-')[1]))
                    self.number_for_label[label] = ('Section', number)
        for page in pages:
            self.page_level = page['level']
            title = self.inline(page['blocks'][0]['c'][2])
            body = self.blocks(page['blocks'])
            (self.dest / (page['name'] + '.md')).write_text(body + '\n')
            self.record['sections'].append({'page': 'tutorial/' + page['name'], 'title': title,
                                            'level': page['level'], 'section': page['section'],
                                            'sha256': sha(self.dest / (page['name'] + '.md'))})
        if self.figure_count != len(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{', source)):
            raise ValueError('Figure count mismatch')
        expected_tables = len(re.findall(r'\\begin\{(?:table|longtable)\}', source))
        if self.table_count != expected_tables:
            raise ValueError(f'Table count mismatch: {self.table_count} != {expected_tables}')
        if [x['text'] for x in self.record['commands']] != self.code:
            raise ValueError('Command blocks differ from source listings')
        self.record['labels'] = sorted(self.known_labels)
        (self.root / 'analysis/supplement-v2/conversion-manifest.json').write_text(json.dumps(self.record, indent=2) + '\n')
        print(f'Converted {len(pages)} pages, {self.figure_count} figures, {self.table_count} tables, {len(self.code)} exact code blocks.')
        return self.record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository', type=Path)
    parser.add_argument('--source', type=Path)
    args = parser.parse_args()
    source = args.source or args.repository / 'analysis/supplement-v2/source/supplementary.tex'
    Converter(args.repository, source).run()


if __name__ == '__main__':
    main()

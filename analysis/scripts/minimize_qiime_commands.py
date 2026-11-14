"""Omit explicit defaults using a captured QIIME action signature, not guesses.

Only parameter options are eligible. Inputs, outputs, metadata, non-default
settings, and controls disabling unwanted default analyses are retained.
"""
import json
import re
import shlex


def options(tokens):
    start = None
    result = []
    for i in range(3, len(tokens)):
        if tokens[i].startswith('--'):
            if start is not None:
                result.append(tokens[start:i])
            start = i
    if start is not None:
        result.append(tokens[start:])
    return result


def prefix(tokens):
    return tokens[:next((i for i in range(3, len(tokens)) if tokens[i].startswith('--')), len(tokens))]


def parameter(option, signature):
    flag = option[0]
    if not flag.startswith('--p-'):
        return None, None
    key = flag[4:].replace('-', '_')
    negative = False
    if key not in signature and key.startswith('no_'):
        key = key[3:]
        negative = True
    if key not in signature:
        raise ValueError('Unknown parameter flag: ' + flag)
    spec = signature[key]
    if spec['type'] == 'Bool' or ('Bool' in spec['type'] and len(option) == 1):
        if len(option) != 1:
            raise ValueError('Unexpected explicit Boolean value: ' + str(option))
        return key, not negative
    if negative:
        raise ValueError('Negated non-Boolean flag: ' + flag)
    values = option[1:]
    if not values:
        raise ValueError('Missing parameter value: ' + flag)
    if spec['type'].startswith(('List[', 'Set[')):
        return key, values
    if len(values) != 1:
        raise ValueError('Unexpected multiple parameter values: ' + str(option))
    value = values[0]
    try:
        if spec['type'].startswith('Int'):
            value = int(value)
        elif spec['type'].startswith('Float'):
            value = float(value)
    except ValueError:
        if not (value.startswith('<') and value.endswith('>')):
            raise
    return key, value


def resolved(tokens, snapshot):
    signature = snapshot['actions'].get('/'.join(tokens[1:3]))
    if signature is None:
        return {'unregistered_cli': tokens}
    result = {key: spec['default'] for key, spec in signature.items() if not spec['required']}
    nonparams = []
    for option in options(tokens):
        key, value = parameter(option, signature)
        if key is None:
            nonparams.append(option)
        else:
            result[key] = value
    return {'parameters': result, 'inputs_outputs_metadata': nonparams, 'command_prefix': prefix(tokens)}


def quote(value):
    if value.startswith('<') or value.endswith('>'):
        return value
    if re.fullmatch(r'[A-Za-z0-9_./:@%+=,-]+', value):
        return value
    if ('${' in value or "'" in value) and '"' not in value:
        return '"' + value + '"'
    return shlex.quote(value)


def render(tokens):
    chunks = [' '.join(quote(t) for t in prefix(tokens))]
    for option in options(tokens):
        line = '    '
        for index, token in enumerate(option):
            word = quote(token)
            if index > 1 and len(line) + len(word) > 88 and line.strip():
                chunks.append(line.rstrip())
                line = '        '
            line += word + ' '
        chunks.append(line.rstrip())
    return ' \\\n'.join(chunks)


def simplify(command, snapshot):
    tokens = shlex.split(command, comments=True)
    if tokens[:1] != ['qiime']:
        return command, None
    signature = snapshot['actions'].get('/'.join(tokens[1:3]))
    if signature is None:
        return command, None
    output = prefix(tokens)
    removed = []
    for option in options(tokens):
        key, value = parameter(option, signature)
        spec = signature.get(key, {})
        if key is not None and not spec.get('required', True) and value == spec['default']:
            removed.append({'flag': option[0], 'parameter': key, 'value': value, 'default': spec['default']})
        else:
            output.extend(option)
    assert resolved(tokens, snapshot) == resolved(output, snapshot), 'Effective parameter settings changed'
    # Reparse the rendered shell, checking quotation and line continuation too.
    new = render(output)
    assert shlex.split(new.replace('\\\n', ' '), comments=True) == output
    return new, {'action': '/'.join(tokens[1:3]), 'before': tokens, 'after': output, 'removed': removed,
                 'effective_settings_equal': True}


def simplify_block(block, snapshot):
    lines = block.splitlines()
    output, records = [], []
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.lstrip().startswith('qiime '):
            output.append(line)
            i += 1
            continue
        parts = [line]
        while parts[-1].rstrip().endswith('\\'):
            i += 1
            if i >= len(lines):
                raise ValueError('Unterminated shell continuation')
            parts.append(lines[i])
        command = '\n'.join(parts).replace('\\\n', ' ')
        replacement, record = simplify(command, snapshot)
        output.append(replacement)
        if record:
            records.append(record)
        i += 1
    return '\n'.join(output), records


def simplify_text(text, snapshot, latex=False):
    records = []
    pattern = (r'(\\begin\{(?:PseudoCode|CodeInput|CodeOutput)\}\n)(.*?)(\n\\end\{(?:PseudoCode|CodeInput|CodeOutput)\})'
               if latex else r'(^```bash\n)(.*?)(\n```$)')
    def replace(match):
        body, found = simplify_block(match.group(2), snapshot)
        records.extend(found)
        return match.group(1) + body + match.group(3)
    return re.sub(pattern, replace, text, flags=re.M | re.S), records

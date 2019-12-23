"""
Functions helping to nicely format values such as farads, ohms, etc.
"""

# scales with their prefixes such as kilo- mega- micro- etc.
from typing import Sequence, Tuple, List
from IPython.display import Markdown

SCALES = ((1e6, 'Meg'), (1e3, 'k'), (1, ''), (1e-3, 'm'), (1e-6, 'u'),
          (1e-9, 'n'), (1e-12, 'p'), (1e-15, 'f'), (1e-18, 'a'))


def value_scale(value) -> (float, str):
    """
    Determines the scale and the proper prefix of a value.
    """
    if value == 0:
        return 1.0, ''
    for scale, prefix in SCALES:
        if abs(value) >= scale:
            return scale, prefix
    return 1.0, 'no prefix'


def format_value(value: float, unit: str) -> str:
    """
    Formats a value using a unit. For example, it turns 1e-5 and 'F' into
    '10 uF'.
    """
    scale, prefix = value_scale(value)
    return str(round(value/scale, 2)) + prefix + unit


# noinspection PyPep8Naming
def format_C(value):
    return format_value(value, 'F')


# noinspection PyPep8Naming
def format_R(value):
    return format_value(value, 'Ω')


# noinspection PyPep8Naming
def format_L(value):
    return format_value(value, 'H')


# noinspection PyPep8Naming
def format_I(value):
    return format_value(value, 'A')


# noinspection PyPep8Naming
def format_V(value):
    return format_value(value, 'V')


# noinspection PyPep8Naming
def format_F(value):
    return format_value(value, 'Hz')


# noinspection PyPep8Naming
def format_W(value):
    return format_value(value, 'W')


def pad_str(s: str, length: int, spacer: str) -> str:
    """
    Adds spaces to the end of a string to make its length equal to `length`.
    To make it work correctly, the spacer must be just 1 character long.
    """
    if len(s) < length:
        return s + spacer * (length - len(s))
    return s


def format_markdown_table(headers: Sequence[str],
                          rows: Sequence[Sequence[str]]) -> str:
    """
    Generates a markdown-style table out of a list of headers and rows.
    """
    col_width = [len(h) for h in headers]
    for row in rows:
        assert len(row) == len(headers)
        for i, c in enumerate(row):
            if len(c) > col_width[i]:
                col_width[i] = len(c)
    header_str = (
            '| '
            + ' | '.join(pad_str(h, col_width[i], ' ')
                         for i, h in enumerate(headers))
            + ' |')
    header_separator = (
            '|-'
            + '-|-'.join(pad_str('-', col_width[i], '-')
                         for i, h in enumerate(headers))
            + '-|')
    table = '\n'.join(
        [header_str, header_separator] +
        [('| '
          + ' | '.join(pad_str(h, col_width[i], ' ')
                       for i, h in enumerate(row))
          + ' |')
         for row in rows])
    return table


def warning_message(text):
    return Markdown(
        f'<span style="font-weight:bold; color: #ff0000;">WARNING:</span> '
        f'<span style="font-weight:bold;">{text}</span>')


PrintableValue = Tuple[str, str]


# a bunch of values that should be printed where each value consists of
# (<description>, <the value itself as a string>)
PrintableValues = List[PrintableValue]


def block_of_values_plain(*data: PrintableValue) -> str:
    output = [
        (f'<tr><td>{note}</td><td>{value}</td></tr>' if note
         else f'<tr><td colspan="2">{value}</td></tr>')
        for note, value in data]
    return '<table>' + '\n'.join(output) + '</table>'


def block_of_values(*data: PrintableValue) -> Markdown:
    return Markdown(block_of_values_plain(*data))


def format_flux_density(value: float, saturation_level: float):
    """
    Formats magnetic flux density and adds a warning if it's too high
    to saturate the core.
    """
    base_message = format_value(value, "T")
    if value >= saturation_level:
        result = (
            base_message
            + ' '
            + warning_message(
                f'Such flux density is too high (>= '
                f'{format_value(saturation_level, "T")}) '
                f'and will saturate the core!').data)
    else:
        result = base_message + ', which is OK'
    return result


def test_format_markdown_table():
    assert format_markdown_table(['h1', 'h2'], [['a', 'b'], ['c', 'd']]) == (
        '| h1 | h2 |\n'
        '|----|----|\n'
        '| a  | b  |\n'
        '| c  | d  |')


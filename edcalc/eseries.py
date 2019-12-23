"""
This module is all about preferred values of electronic components.
For example, it contains functions helping to find a standard resistor value
which is closes to a given one, or calculate a voltage divider which will
give the most accurate output voltage with just two standard resistors,
and so on.
"""

from typing import Sequence, List, Tuple
import itertools
import math

import numpy as np

from .format import value_scale, format_C, format_R

# Basic values of E12 series (will be expanded to include x10, x100, and x1000)
e12_series = (1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2)
# Basic values of E24 series (will be expanded to include x10, x100, and x1000)
e24_series = (1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0,
              3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1)


def _expand_exx_series(series: Sequence):
    """
    Expand series by adding, x10, x100 and x1000 values
    """
    return tuple(
        [series[-1] / 10]
        + list(series)
        + [round(s * 10) for s in series]
        + [round(s * 100) for s in series]
        + [series[0] * 1000])


e12_full_series = _expand_exx_series(e12_series)
e24_full_series = _expand_exx_series(e24_series)


def find_voltage_divider(
        v_in: float, v_out: float,
        r1_series: Sequence[float] = e12_full_series,
        r2_series: Sequence[float] = e12_full_series,
        multiplier: float = 1,
        top: int = 5) -> List[Tuple[float, float]]:
    """
    Finds the best voltage divider made of just two resistors of given
    series (E12/E24).
    :param v_in: Input voltage
    :param v_out: Output voltage
    :param r1_series: Enn series to choose the upper resistor from.
    :param r2_series: Enn series to choose for the lower resistor from.
    :param multiplier: multiply the result by a given factor.
    :param top: how many of the combinations to return
    :return: The result is a list of combinations
    (<R1 upper>, <R2 lower>) with the best one on the top.
    """
    best_ratio = v_in / v_out - 1
    all_combinations = sorted(
        [(abs(best_ratio - r1 / r2) / best_ratio, r1, r2)
         for r1, r2 in itertools.product(r1_series, r2_series)],
        key=lambda x: x[0])
    unique_ratios = all_combinations[:1] + [
        all_combinations[i]
        for i in range(1, len(all_combinations))
        if abs(all_combinations[i][0] - all_combinations[i - 1][0]) > 1e-6]
    result = [(x[0], x[1] * multiplier, x[2] * multiplier)
              for x in unique_ratios]
    return [tuple(x[1:]) for x in result[:top]]


def find_esum(value: float, max_components: int = 2,
              std_values: Sequence = e12_full_series,
              max_difference: float = 1e6) -> List[float]:
    """
    Finds a collection of standard values that can be summed up to get
    the given value, or something close to it.

    :param value: the value we are trying to split
    :param max_components: specifies how many elements the sum should contain.
    :param std_values: a complete series of standard values, like
        `e24_full_series`.
    :param max_difference: order of magnitude of difference between
        the adjacent components of the sum. This means that if we split number
        10000001 with max_difference=1e6, the result will be a single
        value of 1Meg, and not 1Meg + 1. Given 10000011 though, the result
        would be 1Meg + 10 + 1.
    """
    if max_components < 2:
        return [nearest_standard(value, std_values)]
    else:
        scale, _ = value_scale(value)
        standard = nearest_standard(value, std_values)
        if math.isclose(standard, value, abs_tol=scale / max_difference):
            return [standard]
        else:
            lower_standard = nearest_standard(value, std_values, 'lower')
            return ([lower_standard]
                    + find_esum(value - lower_standard,
                                max_components - 1,
                                std_values))


def nearest_standard(value: float, std_values: Sequence,
                     strategy: str = 'any') -> float:
    """
    Finds a closest standard value to a given preferred series like E12 or E24.
    Normally this function is not supposed to be called directly. Use
    `nearest_e12` or `nearest_e24` instead.

    :param value: a precise value we need to replace with a standard one
    :param std_values: a complete series of standard values, like
        `e24_full_series`.
    :param Literal["any", "higher", "lower"] strategy: a rounding strategy,
      where "any" means we're looking for any closest standard value,
      "higher" means the standard value must be larger than the precise one,
      and "lower" means the opposite.
    """
    if value == 0:
        return value
    std = np.array(std_values)
    scale, prefix = value_scale(value)
    dist_2 = (std - value / scale)**2
    match_index = np.argmin(dist_2)
    result = std[match_index] * scale
    if strategy == 'any':
        return result
    elif strategy == 'higher':
        if result >= value:
            return result
        else:
            return std[match_index + 1] * scale
    elif strategy == 'lower':
        if result <= value:
            return result
        else:
            return std[match_index - 1] * scale
    else:
        raise ValueError(f'Unknown rounding strategy: {strategy!r}')


def nearest_e24(value: float, strategy: str = 'any') -> float:
    """
    Finds a value from E24 series which is closest to the given one, with
    some rounding strategy in mind.

    :param value: a precise value we need to replace with a standard one
    :param Literal["any", "higher", "lower"] strategy: a rounding strategy,
      where "any" means we're looking for any closest standard value,
      "higher" means the standard value must be larger than the precise one,
      and "lower" means the opposite.
    """
    return nearest_standard(
        value,
        # the first and the last values are necessary to make all strategies
        # work properly
        e24_full_series,
        strategy)


def nearest_e12(value, strategy: str = 'any'):
    """
    Finds a value from E12 series which is closest to the given one, with
    some rounding strategy in mind.

    :param value: a precise value we need to replace with a standard one
    :param Literal["any", "higher", "lower"] strategy: a rounding strategy,
      where "any" means we're looking for any closest standard value,
      "higher" means the standard value must be larger than the precise one,
      and "lower" means the opposite.
    """
    return nearest_standard(
        value,
        # the first and the last values are necessary to make all strategies
        # work properly
        e12_full_series,
        strategy)


def test_nearest_e24():
    from pytest import approx
    assert nearest_e24(1e-9) == approx(1e-9)
    assert nearest_e24(9.2e-9) == approx(9.1e-9)
    assert nearest_e12(9.2e-9) == approx(1e-8)


def test_find_esum():
    assert [format_C(c) for c in find_esum(11e-6, 2)] == ['10.0uF', '1.0uF']
    assert [format_R(c) for c in find_esum(1000001, 2)] == ['1.0MegΩ']
    assert ([format_R(c) for c in find_esum(1000011, 3)]
            == ['1.0MegΩ', '10.0Ω', '1.0Ω'])

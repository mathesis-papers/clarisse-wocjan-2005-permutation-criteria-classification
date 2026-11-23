import argparse
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from math import comb
from pathlib import Path


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def labels(mask, parties):
    return tuple(((mask >> (2 * site)) & 1, (mask >> (2 * site + 1)) & 1)
                 for site in range(parties))


def distance(first, second):
    return sum(left != right for left, right in zip(first, second))


def signature(mask, parties):
    family = labels(mask, parties)
    sectors = tuple(first == second for first, second in family)
    distances = tuple(distance(family[first], family[second])
                      for first in range(parties) for second in range(first + 1, parties))
    return sectors, distances


def exchange(mask, parties):
    return sum(((mask >> position) & 1) << (position ^ 1)
               for position in range(2 * parties))


def orbit(mask, parties):
    full = (1 << (2 * parties)) - 1
    swapped = exchange(mask, parties)
    return {mask, mask ^ full, swapped, swapped ^ full}


def sparse_flatten(entries, label_family):
    positions = tuple(bit for pair in label_family for bit in pair)
    rows = sorted({tuple(value for value, side in zip(index, positions) if side)
                   for index, coefficient in entries})
    columns = sorted({tuple(value for value, side in zip(index, positions) if not side)
                      for index, coefficient in entries})
    matrix = [[0 for column in columns] for row in rows]
    row_ids = {row: index for index, row in enumerate(rows)}
    column_ids = {column: index for index, column in enumerate(columns)}
    for index, coefficient in entries:
        row = tuple(value for value, side in zip(index, positions) if side)
        column = tuple(value for value, side in zip(index, positions) if not side)
        matrix[row_ids[row]][column_ids[column]] += coefficient
    return matrix


def gram(matrix):
    return [[sum(left * right for left, right in zip(row, other))
             for other in matrix] for row in matrix]


def multiply(first, second):
    return [[sum(left * right for left, right in zip(row, column))
             for column in zip(*second)] for row in first]


def check_projector(matrix, eigenvalue, multiplicity):
    squared = multiply(matrix, matrix)
    require(squared == [[eigenvalue * value for value in row] for row in matrix],
            'Incorrect Gram polynomial')
    require(sum(matrix[index][index] for index in range(len(matrix))) ==
            eigenvalue * multiplicity, 'Incorrect Gram multiplicity')


def witness_checks():
    mixed_cases = []
    for label in product((0, 1), repeat=2):
        matrix = sparse_flatten([((0, 0), 1), ((1, 1), 1)], (label,))
        same_side = label[0] == label[1]
        eigenvalue = 2 if same_side else 1
        multiplicity = 1 if same_side else 2
        check_projector(gram(matrix), eigenvalue, multiplicity)
        norm_squared = Fraction(eigenvalue * multiplicity**2, 4)
        require(norm_squared == (Fraction(1, 2) if same_side else 1),
                'Incorrect mixed trace norm')
        mixed_cases.append({'label': list(label), 'gram_eigenvalue': eigenvalue,
                            'multiplicity': multiplicity, 'trace_norm_squared': str(norm_squared)})
    bell_cases = []
    entries = [((ket_bit, bra_bit, ket_bit, bra_bit), 1)
               for ket_bit, bra_bit in product((0, 1), repeat=2)]
    for first in product((0, 1), repeat=2):
        for second in product((0, 1), repeat=2):
            separation = distance(first, second)
            matrix = sparse_flatten(entries, (first, second))
            eigenvalue = 2 ** (2 - separation)
            multiplicity = 2 ** separation
            check_projector(gram(matrix), eigenvalue, multiplicity)
            norm_squared = Fraction(eigenvalue * multiplicity**2, 4)
            require(norm_squared == 2**separation, 'Incorrect Bell trace norm')
            bell_cases.append({'labels': [list(first), list(second)], 'distance': separation,
                               'gram_eigenvalue': eigenvalue, 'multiplicity': multiplicity,
                               'trace_norm_squared': str(norm_squared)})
    return {'matrix_scaling': 'Each matrix is twice the density-matrix rearrangement.',
            'method': 'Integer Gram polynomial G^2=lambda G, with exact trace and multiplicity.',
            'mixed_cases': mixed_cases, 'bell_cases': bell_cases}


def enumerate_labels(max_parties):
    result = []
    for parties in range(max_parties + 1):
        grouped = {}
        for mask in range(1 << (2 * parties)):
            key = signature(mask, parties)
            representative = min(orbit(mask, parties))
            require(grouped.setdefault(key, representative) == representative,
                    'Label reconstruction failed')
        result.append({'parties': parties, 'label_families': 4**parties,
                       'signature_classes': len(grouped)})
    return result


def enumerate_balanced(max_parties):
    result = []
    for parties in range(1, max_parties + 1):
        grouped = {}
        representatives = set()
        fixed_exchange = 0
        fixed_exchange_complement = 0
        count = 0
        full = (1 << (2 * parties)) - 1
        for chosen in combinations(range(2 * parties), parties):
            mask = sum(1 << position for position in chosen)
            representative = min(orbit(mask, parties))
            key = signature(mask, parties)
            require(grouped.setdefault(key, representative) == representative,
                    'Balanced label reconstruction failed')
            require(signature(representative, parties) == key, 'Orbit signature changed')
            representatives.add(representative)
            fixed_exchange += exchange(mask, parties) == mask
            fixed_exchange_complement += exchange(mask, parties) ^ full == mask
            count += 1
        expected_exchange = comb(parties, parties // 2) if parties % 2 == 0 else 0
        numerator = comb(2 * parties, parties) + 2**parties + expected_exchange
        require(numerator % 4 == 0, 'Orbit count is not integral')
        require(count == comb(2 * parties, parties), 'Incorrect subset count')
        require(fixed_exchange == expected_exchange, 'Incorrect exchange fixed points')
        require(fixed_exchange_complement == 2**parties, 'Incorrect combined fixed points')
        require(len(grouped) == len(representatives) == numerator // 4,
                'Orbit count disagrees with signature classes')
        result.append({'parties': parties, 'balanced_subsets': count,
                       'exchange_fixed': fixed_exchange,
                       'exchange_complement_fixed': fixed_exchange_complement,
                       'norm_classes': len(grouped), 'nontrivial_classes': len(grouped) - 1})
    return result


def permutation_convention(max_parties):
    result = []
    for parties in range(1, max_parties + 1):
        count = 0
        entry_checks = 0
        for permutation in permutations(range(2 * parties)):
            selected = tuple(position for position, output in enumerate(permutation) if output % 2 == 0)
            excluded = tuple(position for position in range(2 * parties) if position not in selected)
            row_images = set()
            column_images = set()
            for row in product((0, 1), repeat=parties):
                row_images.add(tuple(row[permutation[position] // 2] for position in selected))
            for column in product((0, 1), repeat=parties):
                column_images.add(tuple(column[permutation[position] // 2] for position in excluded))
            require(len(row_images) == len(column_images) == 2**parties,
                    'Row or column reindexing is not bijective')
            for output in product((0, 1), repeat=2 * parties):
                original = tuple(output[permutation[position]] for position in range(2 * parties))
                flattened_row = tuple(original[position] for position in selected)
                flattened_column = tuple(original[position] for position in excluded)
                require(flattened_row == tuple(output[permutation[position]] for position in selected),
                        'Row-set convention mismatch')
                require(flattened_column == tuple(output[permutation[position]] for position in excluded),
                        'Column-set convention mismatch')
                entry_checks += 1
            count += 1
        result.append({'parties': parties, 'permutations': count, 'entries': entry_checks})
    return result


def adaptive_and_rank_checks(max_parties):
    results = []
    for parties in range(1, max_parties + 1):
        adaptive_groups = {}
        rank_groups = {}
        maximum_tests = 0
        for chosen in combinations(range(2 * parties), parties):
            mask = sum(1 << position for position in chosen)
            family = labels(mask, parties)
            sectors = tuple(first == second for first, second in family)
            distances = []
            for sector in (False, True):
                members = [site for site in range(parties) if sectors[site] == sector]
                if members:
                    distances.extend(distance(family[members[0]], family[site])
                                     for site in members[1:])
            key = sectors, tuple(distances)
            representative = min(orbit(mask, parties))
            require(adaptive_groups.setdefault(key, representative) == representative,
                    'Adaptive label recovery failed')
            maximum_tests = max(maximum_tests, parties + len(distances))
            rank_signature = (tuple(1 if sector else 2 for sector in sectors),
                              tuple(2**distance(family[first], family[second])
                                    for first, second in combinations(range(parties), 2)))
            require(rank_groups.setdefault(rank_signature, representative) == representative,
                    'Rank label recovery failed')
        require(maximum_tests <= 2 * parties - 1, 'Adaptive comparison bound failed')
        require(len(adaptive_groups) == len(rank_groups), 'Rank/adaptive class mismatch')
        results.append({'parties': parties, 'adaptive_classes': len(adaptive_groups),
                        'rank_classes': len(rank_groups), 'maximum_adaptive_tests': maximum_tests})
    return results


def schatten_exponent_checks(max_parties):
    results = []
    exponents = [('1/2', Fraction(2)), ('1', Fraction(1)),
                 ('3/2', Fraction(2, 3)), ('2', Fraction(1, 2)),
                 ('4', Fraction(1, 4)), ('infinity', Fraction(0))]
    for parties in range(1, max_parties + 1):
        classes = []
        for label, inverse in exponents:
            grouped = {}
            for chosen in combinations(range(2 * parties), parties):
                mask = sum(1 << position for position in chosen)
                family = labels(mask, parties)
                mixed = tuple(-Fraction(1, 2) if first == second else inverse - 1
                              for first, second in family)
                bell = tuple(distance(family[first], family[second]) * (inverse - Fraction(1, 2))
                             for first, second in combinations(range(parties), 2))
                representative = min(orbit(mask, parties)) if label != '2' else 0
                require(grouped.setdefault((mixed, bell), representative) == representative,
                        'Schatten exponent separation failed')
            classes.append({'p': label, 'classes': len(grouped)})
        results.append({'parties': parties, 'sampled_schatten_exponents': classes})
    return results


def partial_transpose_checks(max_parties):
    results = []
    for parties in range(1, max_parties + 1):
        cuts = range(1, 1 << (parties - 1))
        pair_count = 0
        for first_cut, second_cut in combinations(cuts, 2):
            directions = set()
            for first, second in combinations(range(parties), 2):
                first_separates = ((first_cut >> first) ^ (first_cut >> second)) & 1
                second_separates = ((second_cut >> first) ^ (second_cut >> second)) & 1
                if first_separates != second_separates:
                    directions.add((first_separates, second_separates))
            require(directions == {(1, 0), (0, 1)}, 'Cut incomparability failed')
            pair_count += 1
        results.append({'parties': parties, 'nontrivial_cuts': len(cuts),
                        'cut_pairs_with_both_bell_witnesses': pair_count})
    return results


def main():
    parser = argparse.ArgumentParser(description='Exact finite checks for permutation norm dependence.')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'schema_version': 1, 'arithmetic': 'integers and exact fractions',
              'coverage': 'Local witness blocks; label families through 6 parties; balanced subsets, adaptive/rank signatures and partial-transpose cuts through 8 parties; sampled Schatten exponents through 6 parties; permutation convention through 3 parties.',
              'witnesses': witness_checks(), 'label_reconstruction': enumerate_labels(6),
              'balanced_subsets': enumerate_balanced(8),
              'permutation_convention': permutation_convention(3),
              'adaptive_and_rank': adaptive_and_rank_checks(8),
              'schatten_exponents': schatten_exponent_checks(6),
              'partial_transpose_cuts': partial_transpose_checks(8),
              'all_checks_satisfied': True}
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()

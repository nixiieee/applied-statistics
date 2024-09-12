# test_nist.py
import pytest
from scipy.special import erfc as erfc
from math import sqrt as sqrt
from math import fabs as fabs
from scipy.special import gammaincc as gammaincc
from math import floor as floor
from math import log2

with open("binary_data.txt", "r") as f:
    binary_data = f.read().strip()

def approximate_entropy_test(binary_data, m=2):
    def phi(m):
        n = len(binary_data)
        counts = {}
        for i in range(n - m + 1):
            sub_seq = binary_data[i:i + m]
            counts[sub_seq] = counts.get(sub_seq, 0) + 1
        return sum(count / n * log2(count / n) for count in counts.values())

    n = len(binary_data)
    phi_m = phi(m)
    phi_m_plus_1 = phi(m + 1)
    ap_en = phi_m - phi_m_plus_1
    chi_square = 2 * n * (log2(2) - ap_en)
    p_value = gammaincc(2 ** (m - 1), chi_square / 2)
    return p_value >= 0.01

def test_approximate_entropy():
    assert approximate_entropy_test(binary_data), "Approximate Entropy test failed"

def frequency_test(binary_data):
    count = 0
    for bit in binary_data:
        if bit == '0':
            count -= 1
        elif bit == '1':
            count += 1
    sObs = count / sqrt(len(binary_data))
    p_value = erfc(fabs(sObs) / sqrt(2))
    return p_value >= 0.01

def test_frequency():
    assert frequency_test(binary_data), "Frequency test failed"

def runs_test(binary_data):
    runs = 1
    for i in range(1, len(binary_data)):
        if binary_data[i] != binary_data[i - 1]:
            runs += 1
    total_bits = len(binary_data)
    expected_runs = (2 * binary_data.count('1') * binary_data.count('0')) / total_bits + 1
    return abs(runs - expected_runs) / expected_runs >= 0.01

def test_runs():
    assert runs_test(binary_data), "Runs test failed"

def block_frequency_test(binary_data, block_size=128):
    length_of_bit_string = len(binary_data)


    if length_of_bit_string < block_size:
        block_size = length_of_bit_string

    # Compute the number of blocks based on the input given.
    number_of_blocks = floor(length_of_bit_string / block_size)

    if number_of_blocks == 1:
        # For block size M=1, this test degenerates to test 1, the Frequency (Monobit) test.
        return frequency_test(binary_data[0:block_size])

    # Initialized variables
    block_start = 0
    block_end = block_size
    proportion_sum = 0.0

    # Create a for loop to process each block
    for counter in range(number_of_blocks):
        # Partition the input sequence and get the data for block
        block_data = binary_data[block_start:block_end]

        one_count = 0
        for bit in block_data:
            if bit == '1':
                one_count += 1
        pi = one_count / block_size

        proportion_sum += pow(pi - 0.5, 2.0)

        block_start += block_size
        block_end += block_size

    result = 4.0 * block_size * proportion_sum

    p_value = gammaincc(number_of_blocks / 2, result / 2)

def test_block_frequency():
    assert block_frequency_test(binary_data), "Block Frequency test failed"

def longest_run_of_ones_test(binary_data):
    # Определяем длину самой длинной последовательности единиц
    max_run = 0
    current_run = 0
    for bit in binary_data:
        if bit == '1':
            current_run += 1
            if current_run > max_run:
                max_run = current_run
        else:
            current_run = 0
    
    n = len(binary_data)
    
    # Ожидаемая максимальная длина для случайной строки длиной n
    expected_max_run = floor(log2(n))
    
    # Рассчитываем p-value, сравнивая найденный максимальный пробег с ожидаемым
    p_value = gammaincc(expected_max_run, max_run / 2)
    
    return p_value >= 0.01

def test_longest_run_of_ones():
    assert longest_run_of_ones_test(binary_data), "Longest Run of Ones test failed"

def cumulative_sums_test(binary_data):
    cumulative_sum = 0
    max_deviation = 0
    for i, bit in enumerate(binary_data):
        if bit == '1':
            cumulative_sum += 1
        else:
            cumulative_sum -= 1
        max_deviation = max(max_deviation, abs(cumulative_sum))
    
    p_value = gammaincc(len(binary_data) / 2, max_deviation / 2)
    return p_value >= 0.01

def test_cumulative_sums():
    assert cumulative_sums_test(binary_data), "Cumulative Sums test failed"

def serial_test(binary_data, m=2):
    n = len(binary_data)
    
    def count_patterns(pattern_length):
        counts = {}
        for i in range(n - pattern_length + 1):
            sub_seq = binary_data[i:i + pattern_length]
            counts[sub_seq] = counts.get(sub_seq, 0) + 1
        return counts

    phi_m = sum((count / n) ** 2 for count in count_patterns(m).values())
    phi_m_plus_1 = sum((count / n) ** 2 for count in count_patterns(m + 1).values())
    
    delta_phi = phi_m - phi_m_plus_1
    p_value = gammaincc(2 ** (m - 1), delta_phi / 2)
    
    return p_value >= 0.01

def test_serial():
    assert serial_test(binary_data), "Serial test failed"

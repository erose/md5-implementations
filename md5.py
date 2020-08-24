from typing import *
import sys
import math

Bit = int
Byte = int
Word = int

def md5(input_bytes: bytes) -> str:
    bits = bytes_to_bits(input_bytes)

    bits = preprocess(bits)
    preprocessed_bytes = bits_to_bytes(bits)
    mixed_around_bytes = mix_em_around(preprocessed_bytes)

    # Output is as a hex string.
    return "".join(format(byte, '0x') for byte in mixed_around_bytes)

def reinterpret(data: List[int], old_width: int, new_width: int):
    bits : List[Bit] = []
    for datum in data:
        for i in reversed(range(old_width)):
            mask = 1 << i
            bits.append((mask & datum) >> i)

    if len(bits) % new_width != 0:
        raise ValueError(f"len(bits) must be a multiple of {new_width}, but it is {len(bits)}")

    result = []
    for i in range(0, len(bits), new_width):
        new_datum = sum(bits[i+j] << (new_width - 1 - j) for j in range(new_width))
        result.append(new_datum)

    return result

bits_to_bytes = lambda bits: reinterpret(bits, 1, 8)
bits_to_word = lambda bits: reinterpret(bits, 1, 32)[0]
bytes_to_bits = lambda bytes_: reinterpret(bytes_, 8, 1)
word_to_bits = lambda word: reinterpret([word], 32, 1)

def bytes_to_words(data: List[int]):
   # We need to switch endianness here because "a sequence of
   # bytes can be interpreted as a sequence of 32-bit words, where each
   # consecutive group of four bytes is interpreted as a word with the
   # low-order (least significant) byte given first."
    bytes_in_low_endian_order = sum([list(reversed(data[i:i+4])) for i in range(0, len(data), 4)], [])
    return reinterpret(bytes_in_low_endian_order, 8, 32)

def words_to_bytes(data: List[int]):
    bytes_in_low_endian_order = reinterpret(data, 32, 8)
    # We need to switch endianness here because "a sequence of
    # bytes can be interpreted as a sequence of 32-bit words, where each
    # consecutive group of four bytes is interpreted as a word with the
    # low-order (least significant) byte given first."
    bytes_ = sum([list(reversed(bytes_in_low_endian_order[i:i+4])) for i in range(0, len(bytes_in_low_endian_order), 4)], [])
    return bytes_

def rotate_word_left(x: Word, steps: int):
    bits = word_to_bits(x)

    for _ in range(steps):
        bits.append(bits.pop(0))

    return bits_to_word(bits)

"""
The message is "padded" (extended) so that its length (in bits) is
congruent to 448, modulo 512. That is, the message is extended so
that it is just 64 bits shy of being a multiple of 512 bits long.
Padding is always performed, even if the length of the message is
already congruent to 448, modulo 512.

Padding is performed as follows: a single "1" bit is appended to the
message, and then "0" bits are appended so that the length in bits of
the padded message becomes congruent to 448, modulo 512. In all, at
least one bit and at most 512 bits are appended.
"""
def pad(bits: List[Bit]) -> List[Bit]:
    bits_to_add = 448 - (len(bits) % 512)

    bits.append(1)
    bits_to_add -= 1

    bits.extend([0]*bits_to_add)

    return bits

"""
A 64-bit representation of b (the length of the message before the
padding bits were added) is appended to the result of the previous
step. In the unlikely event that b is greater than 2^64, then only
the low-order 64 bits of b are used. (These bits are appended as two
32-bit words and appended low-order word first in accordance with the
previous conventions.)
"""
def length_as_bits(b: int) -> List[Bit]:
    result: List[Bit] = []
    for i in range(64):
        result.append((b & (1 << i)) >> i)

    return result

def preprocess(bits: List[Bit]) -> List[Bit]:
    original_length = len(bits)

    bits = pad(bits)
    bits += length_as_bits(original_length)
    
    return bits

# This step uses a 64-element table T[1 ... 64] constructed from the
# sine function. Let T[i] denote the i-th element of the table, which
# is equal to the integer part of 4294967296 times abs(sin(i)), where i
# is in radians.
def construct_T() -> Dict[int, int]:
    return { i: int(4294967296*abs(math.sin(i))) for i in range(1, 64 + 1) }

def mix_em_around(bytes_: List[Byte]) -> List[Byte]:
    # A four-word buffer (A,B,C,D) is used to compute the message digest.
    # Here each of A, B, C, D is a 32-bit register. These registers are
    # initialized to the following values in hexadecimal, low-order bytes
    # first):
    # 
    #    word A: 01 23 45 67
    #    word B: 89 ab cd ef
    #    word C: fe dc ba 98
    #    word D: 76 54 32 10
    A: Word = 0x67452301
    B: Word = 0xefcdab89
    C: Word = 0x98badcfe
    D: Word = 0x10325476

    # We first define four auxiliary functions that each take as input
    # three 32-bit words and produce as output one 32-bit word.
    F = lambda x, y, z: x&y | (~x)&z
    G = lambda x, y, z: x&z | y&(~z)
    H = lambda x, y, z: x ^ y ^ z
    I = lambda x, y, z: y ^ (x | ~z)

    # This step uses a 64-element table T[1 ... 64] constructed from the
    # sine function.
    T = construct_T()

    # At this point the resulting message (after padding with bits and with
    # b) has a length that is an exact multiple of 512 bits. Equivalently,
    # this message has a length that is an exact multiple of 16 (32-bit)
    # words. Let M[0 ... N-1] denote the words of the resulting message,
    # where N is a multiple of 16.
    words: List[Word] = bytes_to_words(bytes_)

    for i in range(len(words) // 16):
        X: List[Word] = [0]*16

        for j in range(16):
            X[j] = words[i*16 + j]

        AA = A
        BB = B
        CC = C
        DD = D

        # /* Round 1. */
        # /* Let [abcd k s i] denote the operation
        #      a = b + ((a + F(b,c,d) + X[k] + T[i]) <<< s). */
        def round1_op(a, b, c, d, k, s, t):
            a += F(b, c, d) + X[k] + T[t]
            a %= 1 << 32
            a = rotate_word_left(a, s)
            a += b
            a %= 1 << 32

            return a

        A = round1_op(A, B, C, D,  0,  7,  1)
        D = round1_op(D, A, B, C,  1, 12,  2)
        C = round1_op(C, D, A, B,  2, 17,  3)
        B = round1_op(B, C, D, A,  3, 22,  4)
        
        A = round1_op(A, B, C, D,  4,  7,  5)
        D = round1_op(D, A, B, C,  5, 12,  6)
        C = round1_op(C, D, A, B,  6, 17,  7)
        B = round1_op(B, C, D, A,  7, 22,  8)

        A = round1_op(A, B, C, D,  8,  7,  9)
        D = round1_op(D, A, B, C,  9, 12, 10)
        C = round1_op(C, D, A, B, 10, 17, 11)
        B = round1_op(B, C, D, A, 11, 22, 12)
        
        A = round1_op(A, B, C, D, 12,  7, 13)
        D = round1_op(D, A, B, C, 13, 12, 14)
        C = round1_op(C, D, A, B, 14, 17, 15)
        B = round1_op(B, C, D, A, 15, 22, 16)

        # /* Round 2. */
        # /* Let [abcd k s t] denote the operation
        #   a = b + ((a + G(b,c,d) + X[k] + T[t]) <<< s). */
        def round2_op(a, b, c, d, k, s, t):
            a += G(b, c, d) + X[k] + T[t]
            a %= 1 << 32
            a = rotate_word_left(a, s)
            a += b
            a %= 1 << 32

            return a

        # /* Do the following 16 operations. */
        A = round2_op(A, B, C, D,  1,  5, 17)
        D = round2_op(D, A, B, C,  6,  9, 18)
        C = round2_op(C, D, A, B, 11, 14, 19)
        B = round2_op(B, C, D, A,  0, 20, 20)
        
        A = round2_op(A, B, C, D,  5,  5, 21)
        D = round2_op(D, A, B, C, 10,  9, 22)
        C = round2_op(C, D, A, B, 15, 14, 23)
        B = round2_op(B, C, D, A,  4, 20, 24)

        A = round2_op(A, B, C, D,  9,  5, 25)
        D = round2_op(D, A, B, C, 14,  9, 26)
        C = round2_op(C, D, A, B,  3, 14, 27)
        B = round2_op(B, C, D, A,  8, 20, 28)
        
        A = round2_op(A, B, C, D, 13,  5, 29)
        D = round2_op(D, A, B, C,  2,  9, 30)
        C = round2_op(C, D, A, B,  7, 14, 31)
        B = round2_op(B, C, D, A, 12, 20, 32)

        # /* Round 3. */
        # /* Let [abcd k s t] denote the operation
        #      a = b + ((a + H(b,c,d) + X[k] + T[t]) <<< s). */
        def round3_op(a, b, c, d, k, s, t):
            a += H(b, c, d) + X[k] + T[t]
            a %= 1 << 32
            a = rotate_word_left(a, s)
            a += b
            a %= 1 << 32

            return a

        # /* Do the following 16 operations. */
        A = round3_op(A, B, C, D,  5,  4, 33)
        D = round3_op(D, A, B, C,  8, 11, 34)
        C = round3_op(C, D, A, B, 11, 16, 35)
        B = round3_op(B, C, D, A, 14, 23, 36)

        A = round3_op(A, B, C, D,  1,  4, 37)
        D = round3_op(D, A, B, C,  4, 11, 38)
        C = round3_op(C, D, A, B,  7, 16, 39)
        B = round3_op(B, C, D, A, 10, 23, 40)

        A = round3_op(A, B, C, D, 13,  4, 41)
        D = round3_op(D, A, B, C,  0, 11, 42)
        C = round3_op(C, D, A, B,  3, 16, 43)
        B = round3_op(B, C, D, A,  6, 23, 44)

        A = round3_op(A, B, C, D,  9,  4, 45)
        D = round3_op(D, A, B, C, 12, 11, 46)
        C = round3_op(C, D, A, B, 15, 16, 47)
        B = round3_op(B, C, D, A,  2, 23, 48)

        # /* Round 4. */
        # /* Let [abcd k s t] denote the operation
        #      a = b + ((a + I(b,c,d) + X[k] + T[t]) <<< s). */
        def round4_op(a, b, c, d, k, s, t):
            a += I(b, c, d) + X[k] + T[t]
            a %= 1 << 32
            a = rotate_word_left(a, s)
            a += b
            a %= 1 << 32

            return a

        # /* Do the following 16 operations. */
        A = round4_op(A, B, C, D,  0,  6, 49)
        D = round4_op(D, A, B, C,  7, 10, 50)
        C = round4_op(C, D, A, B, 14, 15, 51)
        B = round4_op(B, C, D, A,  5, 21, 52)

        A = round4_op(A, B, C, D, 12,  6, 53)
        D = round4_op(D, A, B, C,  3, 10, 54)
        C = round4_op(C, D, A, B, 10, 15, 55)
        B = round4_op(B, C, D, A,  1, 21, 56)

        A = round4_op(A, B, C, D,  8,  6, 57)
        D = round4_op(D, A, B, C, 15, 10, 58)
        C = round4_op(C, D, A, B,  6, 15, 59)
        B = round4_op(B, C, D, A, 13, 21, 60)

        A = round4_op(A, B, C, D,  4,  6, 61)
        D = round4_op(D, A, B, C, 11, 10, 62)
        C = round4_op(C, D, A, B,  2, 15, 63)
        B = round4_op(B, C, D, A,  9, 21, 64)

        A = (A + AA) % (1 << 32)
        B = (B + BB) % (1 << 32)
        C = (C + CC) % (1 << 32)
        D = (D + DD) % (1 << 32)

    # The message digest produced as output is A, B, C, D. That is, we
    # begin with the low-order byte of A, and end with the high-order byte
    # of D.
    return words_to_bytes([A, B, C, D])

import unittest

class TestMd5(unittest.TestCase):
    def test_from_bits_inverts_to_bits(self):
        data = list(b"foo")
        self.assertEqual(bits_to_bytes(bytes_to_bits(data)), data)
        self.assertEqual(bits_to_word(word_to_bits(127)), 127)

    def test_padding(self):
        initial = [0, 0, 0, 1, 0, 1, 1, 1]
        expected = initial + [1] + [0]*(448 - 9)

        computed = pad(initial)
        self.assertEqual(len(computed), 448)
        self.assertEqual(computed, expected)

    def test_rotate_word_left(self):
        self.assertEqual(rotate_word_left(0b01, 1), 0b10)
        self.assertEqual(rotate_word_left(0x12345678, 4), 0x23456781)
        self.assertEqual(rotate_word_left(0x12345678, 16), 0x56781234)

    def test_construct_T(self):
        T = construct_T()

        self.assertEqual(T[1], 0xd76aa478)
        self.assertEqual(T[2], 0xe8c7b756)
        self.assertEqual(T[64], 0xeb86d391)

    def test_md5(self):
        self.assertEqual(md5(b"foo"), "acbd18db4cc2f85cedef654fccc4a4d8")

if __name__ == "__main__":
    # unittest.main()
    input_ = sys.stdin.buffer.read()
    print(md5(input_))

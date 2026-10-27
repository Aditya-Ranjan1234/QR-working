# polynomial.py
# Polynomial arithmetic over GF(256)

from .gf256 import gf_add, gf_sub, gf_mul, gf_div

def poly_add(p1: list[int], p2: list[int]) -> list[int]:
    """Add two polynomials over GF(256)."""
    # Pad the shorter one with zeros
    r = [0] * max(len(p1), len(p2))
    for i in range(len(p1)):
        r[i + len(r) - len(p1)] ^= p1[i]
    for i in range(len(p2)):
        r[i + len(r) - len(p2)] ^= p2[i]
    return r

def poly_mul(p1: list[int], p2: list[int]) -> list[int]:
    """Multiply two polynomials over GF(256)."""
    r = [0] * (len(p1) + len(p2) - 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            r[i + j] ^= gf_mul(p1[i], p2[j])
    return r

def poly_eval(p: list[int], x: int) -> int:
    """Evaluate polynomial p at x using Horner's method."""
    y = p[0]
    for i in range(1, len(p)):
        y = gf_add(gf_mul(y, x), p[i])
    return y

def poly_div(num: list[int], den: list[int]) -> tuple[list[int], list[int]]:
    """
    Divide polynomial num by den over GF(256).
    Returns (quotient, remainder).
    """
    # Copy the numerator
    r = list(num)
    
    # Loop over coefficients to divide
    for i in range(len(num) - len(den) + 1):
        coef = r[i]
        if coef != 0:
            for j in range(1, len(den)):
                # XOR subtraction is XOR addition
                r[i + j] ^= gf_mul(den[j], coef)
    
    # The remainder is the part that cannot be divided further
    separator = len(num) - len(den) + 1
    return r[:separator], r[separator:]

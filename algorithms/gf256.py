# gf256.py
# Galois Field GF(2^8) arithmetic implementation

# Primitive polynomial: x^8 + x^4 + x^3 + x^2 + 1 = 285 (0x11d)
PRIMITIVE_POLYNOMIAL = 285

# Exponent and Log tables for alpha = 2
gf_exp = [0] * 512
gf_log = [0] * 256

def init_tables():
    """Precompute the log and exp tables for fast multiplication and division."""
    global gf_exp, gf_log
    x = 1
    for i in range(255):
        gf_exp[i] = x
        gf_log[x] = i
        # Multiply by 2 (shift left by 1)
        x <<= 1
        if x & 0x100:
            x ^= PRIMITIVE_POLYNOMIAL
    
    # Duplicate exp table to avoid doing modulo 255 operations during multiplication
    for i in range(255, 512):
        gf_exp[i] = gf_exp[i - 255]

# Initialize the tables on load
init_tables()

def gf_add(a: int, b: int) -> int:
    """Addition in GF(256) is bitwise XOR."""
    return a ^ b

def gf_sub(a: int, b: int) -> int:
    """Subtraction in GF(256) is identical to addition (bitwise XOR)."""
    return a ^ b

def gf_mul(a: int, b: int) -> int:
    """Multiply two GF(256) elements using precomputed log and exp tables."""
    if a == 0 or b == 0:
        return 0
    return gf_exp[gf_log[a] + gf_log[b]]

def gf_inv(a: int) -> int:
    """Find the multiplicative inverse of a. a * gf_inv(a) == 1."""
    if a == 0:
        raise ZeroDivisionError("GF(256) division by zero.")
    return gf_exp[255 - gf_log[a]]

def gf_div(a: int, b: int) -> int:
    """Divide a by b in GF(256)."""
    if b == 0:
        raise ZeroDivisionError("GF(256) division by zero.")
    if a == 0:
        return 0
    return gf_exp[gf_log[a] + 255 - gf_log[b]]

def gf_pow(a: int, power: int) -> int:
    """Raise a to the power of power in GF(256)."""
    if a == 0:
        return 0 if power > 0 else 1 # Simple convention
    return gf_exp[(gf_log[a] * power) % 255]

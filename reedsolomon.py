# reedsolomon.py
# Pure Python Reed-Solomon Encoder and Decoder over GF(256)
# Based on the verified academic algorithms for GF(2^8)

from gf256 import gf_add, gf_sub, gf_mul, gf_div, gf_pow, gf_inv, gf_exp, gf_log

# --- Polynomial helpers specifically for RS ---

def gf_poly_mul(p, q):
    """Multiply two polynomials over GF(256) (lowest power first)."""
    r = [0] * (len(p) + len(q) - 1)
    for j in range(len(q)):
        for i in range(len(p)):
            r[i + j] ^= gf_mul(p[i], q[j])
    return r

def gf_poly_eval(p, x):
    """Evaluate polynomial p at x (lowest power first)."""
    # p[0] + p[1]*x + p[2]*x^2 + ...
    y = 0
    xp = 1
    for coef in p:
        y ^= gf_mul(coef, xp)
        xp = gf_mul(xp, x)
    return y

# --- RS ENCODER ---

def rs_generator_poly(nsym):
    """Generate generator polynomial G(x) = (x - a^0)(x - a^1)...(x - a^(nsym-1))"""
    g = [1]
    for i in range(nsym):
        # Multiply by (x - a^i) => (a^i + x)
        # Lowest power first: [gf_exp[i], 1]
        g = gf_poly_mul(g, [gf_exp[i], 1])
    return g

def rs_encode_msg(msg, nsym):
    """Encode msg using systematic Reed-Solomon code (returns msg + parity)."""
    gen = rs_generator_poly(nsym)
    # Pad message with zeros to make space for parity (systematic code)
    # The data bytes are the high order coefficients, parity is low order
    # Let's do polynomial division
    # Copy message
    msg_out = [0] * nsym + msg[::-1]
    
    # Division loop
    for i in range(len(msg) - 1, -1, -1):
        coef = msg_out[i + nsym]
        if coef != 0:
            for j in range(len(gen) - 1):
                msg_out[i + j] ^= gf_mul(gen[j], coef)
                
    # Extract parity (first nsym elements)
    parity = msg_out[:nsym]
    # Return message + parity (high order coefficients first, i.e., msg then parity)
    # Wait, QR standard: message is transmitted first, then parity bytes.
    # Our msg_out has parity in msg_out[:nsym], msg in msg_out[nsym:]
    # So we return list(msg) + parity[::-1] or parity?
    # Let's return list(msg) + parity[::-1] because we traversed backwards.
    # Let's just check: systematic remainder is msg_out[:nsym].
    return list(msg) + msg_out[:nsym][::-1]

# --- RS DECODER ---

def rs_calc_syndromes(msg, nsym):
    """Calculate syndromes. S_i = msg(a^i) for i = 0..nsym-1."""
    # Since msg is high-power-first, let's reverse it to evaluate lowest-power-first
    msg_rev = msg[::-1]
    synd = [0] * nsym
    for i in range(nsym):
        synd[i] = gf_poly_eval(msg_rev, gf_exp[i])
    return [0] + synd # 1-indexed for convenience

def rs_find_error_locator(synd, nsym, erase_loc=None):
    """
    Berlekamp-Massey algorithm to compute error locator polynomial Lambda(x).
    """
    erase_loc = erase_loc or []
    # Initialize locator Lambda(x) = 1
    lam = [1]
    # If erasures are known, initialize Lambda to the erasure locator polynomial
    for el in erase_loc:
        X = gf_exp[el]
        lam = gf_poly_mul(lam, [X, 1]) # (1 + X*x) in lowest-power-first
        
    B = list(lam)
    L = len(erase_loc)
    
    for r in range(L + 1, nsym + 1):
        # Calculate discrepancy d
        d = synd[r]
        for j in range(1, len(lam)):
            if r - j > 0:
                d ^= gf_mul(lam[j], synd[r - j])
                
        B = [0] + B
        
        if d != 0:
            if 2 * L < r + len(erase_loc):
                old_lam = list(lam)
                # lam = lam + d * B
                scaled_B = [gf_mul(x, d) for x in B]
                # Pad lam to match scaled_B size
                if len(lam) < len(scaled_B):
                    lam += [0] * (len(scaled_B) - len(lam))
                for j in range(len(scaled_B)):
                    lam[j] ^= scaled_B[j]
                L = r + len(erase_loc) - L
                # B = old_lam / d
                inv_d = gf_inv(d)
                B = [gf_mul(x, inv_d) for x in old_lam]
            else:
                scaled_B = [gf_mul(x, d) for x in B]
                if len(lam) < len(scaled_B):
                    lam += [0] * (len(scaled_B) - len(lam))
                for j in range(len(scaled_B)):
                    lam[j] ^= scaled_B[j]
                    
    # Clean leading zeros from the end (since lowest-power-first, trailing zeros represent higher degree terms that are 0)
    while len(lam) > 0 and lam[-1] == 0:
        lam.pop()
        
    return lam

def rs_find_error_locator_verbose(synd, nsym, erase_loc=None):
    """Verbose Berlekamp-Massey tracking registers for step-by-step math trace."""
    erase_loc = erase_loc or []
    lam = [1]
    for el in erase_loc:
        X = gf_exp[el]
        lam = gf_poly_mul(lam, [X, 1])
        
    B = list(lam)
    L = len(erase_loc)
    
    steps = []
    steps.append({
        "r": 0,
        "d": 0,
        "L": L,
        "lam": list(lam),
        "B": list(B),
        "description": "Initial state (after erasures locator polynomial initialized)" if erase_loc else "Initial state"
    })
    
    for r in range(L + 1, nsym + 1):
        d = synd[r]
        for j in range(1, len(lam)):
            if r - j > 0:
                d ^= gf_mul(lam[j], synd[r - j])
                
        B = [0] + B
        desc = f"r={r}: calculated discrepancy d={d}. "
        
        if d != 0:
            if 2 * L < r + len(erase_loc):
                old_lam = list(lam)
                scaled_B = [gf_mul(x, d) for x in B]
                if len(lam) < len(scaled_B):
                    lam += [0] * (len(scaled_B) - len(lam))
                for j in range(len(scaled_B)):
                    lam[j] ^= scaled_B[j]
                L = r + len(erase_loc) - L
                inv_d = gf_inv(d)
                B = [gf_mul(x, inv_d) for x in old_lam]
                desc += f"d != 0 and 2*L < r + erasures. Lambda updated. L updated to {L}. B updated to old_Lambda / d."
            else:
                scaled_B = [gf_mul(x, d) for x in B]
                if len(lam) < len(scaled_B):
                    lam += [0] * (len(scaled_B) - len(lam))
                for j in range(len(scaled_B)):
                    lam[j] ^= scaled_B[j]
                desc += "d != 0 but 2*L >= r + erasures. Lambda updated. L unchanged. B shifted."
        else:
            desc += "d == 0. No update to Lambda. L unchanged. B shifted."
            
        lam_clean = list(lam)
        while len(lam_clean) > 0 and lam_clean[-1] == 0:
            lam_clean.pop()
            
        steps.append({
            "r": r,
            "d": d,
            "L": L,
            "lam": lam_clean,
            "B": list(B),
            "description": desc
        })
        
    while len(lam) > 0 and lam[-1] == 0:
        lam.pop()
        
    return lam, steps

def rs_find_errors(err_loc, msg_len):
    """Chien search to find root positions of err_loc (Lambda)."""
    errors = []
    # Lambda(x) = 1 + l_1*x + l_2*x^2 + ...
    # We evaluate at x = a^-i for i = 0..msg_len-1
    for i in range(msg_len):
        Xi_inv = gf_exp[255 - i]
        val = gf_poly_eval(err_loc, Xi_inv)
        if val == 0:
            # Position from the start is msg_len - 1 - i
            errors.append(msg_len - 1 - i)
    return errors

def rs_forney_magnitudes(synd, err_loc, err_pos, msg_len):
    """Forney algorithm to compute error magnitudes."""
    nsym = len(synd) - 1
    # 1. Compute Omega(x) = S(x) * Lambda(x) mod x^nsym
    # S(x) starts from S_1
    S = synd[1:]
    omega = gf_poly_mul(S, err_loc)
    # Keep only the first nsym coefficients (modulo x^nsym, lowest powers first)
    omega = omega[:nsym]
    
    # 2. Compute formal derivative of Lambda(x)
    # lam_prime(x) = lam_1 + 2*lam_2*x + 3*lam_3*x^2 + ...
    # In GF(2), coefficients with odd indices remain, others are 0, shifted down
    lam_prime = []
    for i in range(1, len(err_loc)):
        if i % 2 == 1:
            lam_prime.append(err_loc[i])
        else:
            lam_prime.append(0)
            
    magnitudes = []
    for pos in err_pos:
        # Xi = a^(msg_len - 1 - pos)
        # Xi_inv = a^-(msg_len - 1 - pos) = a^(255 - (msg_len - 1 - pos))
        Xi = gf_exp[msg_len - 1 - pos]
        Xi_inv = gf_inv(Xi)
        
        # Evaluate Omega(Xi_inv)
        num = gf_poly_eval(omega, Xi_inv)
        # Evaluate Lambda'(Xi_inv)
        denom = gf_poly_eval(lam_prime, Xi_inv)
        
        # Magnitude: Y_i = Omega(Xi_inv) / Lambda'(Xi_inv)
        # Sometimes standard formula needs an extra factor depending on roots
        # For roots starting at a^0 (which we do): Y_i = Omega(Xi_inv) / Lambda'(Xi_inv)
        mag = gf_mul(Xi, gf_div(num, denom))
        magnitudes.append(mag)
        
    return magnitudes

def rs_correct_msg(msg, nsym, erase_loc=None):
    """Full decode and repair pipeline."""
    erase_loc = erase_loc or []
    corrected = list(msg)
    
    # Calculate syndromes
    synd = rs_calc_syndromes(msg, nsym)
    if max(synd) == 0:
        # No errors
        return corrected, [], []
        
    # Translate erasure index from start to index from end for polynomial evaluation
    # pos_from_end = len(msg) - 1 - pos
    # BM uses the position from the end!
    erase_loc_end = [len(msg) - 1 - e for e in erase_loc]
    
    # Find Lambda(x)
    err_loc = rs_find_error_locator(synd, nsym, erase_loc_end)
    
    # Find error positions (indices from the start)
    err_pos = rs_find_errors(err_loc, len(msg))
    
    if not err_pos:
        raise ValueError("Could not locate errors (too many errors).")
        
    if len(err_pos) * 2 - len(erase_loc) > nsym:
        raise ValueError("Too many errors/erasures to correct.")
        
    # Calculate magnitudes
    magnitudes = rs_forney_magnitudes(synd, err_loc, err_pos, len(msg))
    
    # Apply corrections
    for pos, mag in zip(err_pos, magnitudes):
        corrected[pos] ^= mag
        
    # Verify syndromes are now all 0
    final_synd = rs_calc_syndromes(corrected, nsym)
    if max(final_synd) != 0:
        raise ValueError("Decoder failed to repair all errors.")
        
    return corrected, err_pos, magnitudes

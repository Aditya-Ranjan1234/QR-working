# test_qr.py
# Verification script for GF(256), Reed-Solomon and QR Encoder/Decoder

from gf256 import gf_mul, gf_div, gf_inv, gf_add
from reedsolomon import rs_encode_msg, rs_correct_msg
from qr_encoder import create_qr_code
from qr_decoder import decode_qr_grid

def test_gf256():
    print("Testing GF(256)...")
    # Multiplicative Inverses: a * inv(a) == 1
    for a in range(1, 256):
        inv_a = gf_inv(a)
        prod = gf_mul(a, inv_a)
        assert prod == 1, f"Inverse failed for {a}: {a} * {inv_a} = {prod}"
    
    # Division
    assert gf_div(10, 10) == 1
    assert gf_div(100, 5) == gf_mul(100, gf_inv(5))
    print("GF(256) passed.")

def test_reedsolomon():
    print("Testing Reed-Solomon...")
    msg = [65, 66, 67, 68, 69, 70]
    nsym = 6
    encoded = rs_encode_msg(msg, nsym)
    assert len(encoded) == len(msg) + nsym
    
    # Decodes perfectly with no errors
    corrected, err_pos, mags = rs_correct_msg(encoded, nsym)
    assert corrected == encoded
    assert len(err_pos) == 0
    
    # Decodes with errors (max capacity = nsym // 2 = 3 errors)
    corrupted = list(encoded)
    corrupted[1] ^= 45
    corrupted[3] ^= 99
    corrupted[7] ^= 12
    
    corrected, err_pos, mags = rs_correct_msg(corrupted, nsym)
    assert corrected == encoded
    assert sorted(err_pos) == [1, 3, 7]
    print("Reed-Solomon passed.")

def test_qr_pipeline():
    print("Testing QR Pipeline...")
    text = "https://example.com"
    matrix, mask, full_bytes = create_qr_code(text, 'M')
    
    decoded_text, err_pos, mags = decode_qr_grid(matrix.grid, matrix.types)
    assert decoded_text == text, f"Expected {text}, got {decoded_text}"
    print("QR Pipeline passed.")

if __name__ == '__main__':
    test_gf256()
    test_reedsolomon()
    test_qr_pipeline()
    print("All tests passed successfully!")

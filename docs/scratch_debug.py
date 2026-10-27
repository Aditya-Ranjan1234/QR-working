# scratch_debug.py
from reedsolomon import rs_encode_msg, rs_correct_msg, rs_calc_syndromes, rs_find_error_locator, rs_find_errors, rs_forney_magnitudes

msg = [65, 66, 67, 68, 69, 70]
nsym = 6
encoded = rs_encode_msg(msg, nsym)

corrupted = list(encoded)
corrupted[1] ^= 45
corrupted[3] ^= 99
corrupted[7] ^= 12

print("Encoded:  ", encoded)
print("Encoded Syndromes:", rs_calc_syndromes(encoded, nsym))
print("Corrupted:", corrupted)

synd = rs_calc_syndromes(corrupted, nsym)
print("Syndromes:", synd)

err_loc = rs_find_error_locator(synd, nsym)
print("Err Locator Poly Lambda(x):", err_loc)

err_pos = rs_find_errors(err_loc, len(corrupted))
print("Chien Search Locations:", err_pos)

if err_pos:
    magnitudes = rs_forney_magnitudes(synd, err_loc, err_pos, len(corrupted))
    print("Forney Magnitudes:", magnitudes)
    
    corrected = list(corrupted)
    for pos, mag in zip(err_pos, magnitudes):
        corrected[pos] ^= mag
    print("Corrected:", corrected)
    print("Matches Encoded:", corrected == encoded)

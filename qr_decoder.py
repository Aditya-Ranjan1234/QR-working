# qr_decoder.py
# Pure Python QR Code Decoder & Error Correction pipeline (Supports Version 1 & 2)

from gf256 import gf_add
from reedsolomon import rs_correct_msg
from qr_encoder import QR_SPECS, EC_INDICATORS, FORMAT_MASK, get_format_info, bits_to_bytes

def decode_format_info(left_bits: list[int], right_bits: list[int]) -> tuple[str, int]:
    """
    Decode the 15-bit format info by comparing with all 32 possible combinations.
    Uses Hamming distance to choose the most likely correct format info, allowing recovery
    even if format modules themselves are damaged.
    """
    best_ec = 'M'
    best_mask = 0
    best_dist = 999
    
    # Check both left and right copies (or combine them)
    # We'll compare each copy individually to all valid code words and pick the minimum distance.
    for ec in ['L', 'M', 'Q', 'H']:
        for mask in range(8):
            valid_bits = get_format_info(ec, mask)
            
            # Calculate Hamming distance for left copy
            dist_l = sum(1 for b1, b2 in zip(left_bits, valid_bits) if b1 != b2)
            # Calculate Hamming distance for right copy
            dist_r = sum(1 for b1, b2 in zip(right_bits, valid_bits) if b1 != b2)
            
            min_dist = min(dist_l, dist_r)
            if min_dist < best_dist:
                best_dist = min_dist
                best_ec = ec
                best_mask = mask
                
    return best_ec, best_mask

def extract_format_bits(grid: list[list[int]]) -> tuple[list[int], list[int]]:
    """Extract format bits from the two redundant locations on the grid."""
    size = len(grid)
    
    # Left copy coordinates
    left_coords = [
        (0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (7, 8), (8, 8),
        (8, 7), (8, 5), (8, 4), (8, 3), (8, 2), (8, 1), (8, 0)
    ]
    
    # Right/bottom copy coordinates
    right_coords = []
    for i in range(8):
        right_coords.append((8, size - 1 - i))
    for i in range(7):
        right_coords.append((size - 7 + i, 8))
        
    left_bits = [grid[r][c] if grid[r][c] is not None else 0 for r, c in left_coords]
    right_bits = [grid[r][c] if grid[r][c] is not None else 0 for r, c in right_coords]
    
    return left_bits, right_bits

def extract_and_unmask_data_bits(grid: list[list[int]], types: list[list[str]], mask_pattern: int) -> tuple[list[int], list[int]]:
    """
    Traverse the grid in a zig-zag path, unmasking only the data modules.
    If grid[r][c] is -1 (scratched), it marks that codeword byte as an erasure.
    Returns (bitstream, erasure_bytes).
    """
    size = len(grid)
    bitstream = []
    erasure_bits = []
    
    col = size - 1
    upward = True
    bit_idx = 0
    
    while col > 0:
        if col == 6:
            col -= 1
            
        rows = range(size - 1, -1, -1) if upward else range(size)
        for r in rows:
            for c in (col, col - 1):
                if types[r][c] == 'data':
                    val = grid[r][c]
                    
                    if val == -1:
                        # Scratched module
                        erasure_bits.append(bit_idx)
                        val = 0 # treat as 0 for bit extraction
                    elif val is None:
                        val = 0
                        
                    # Apply unmasking
                    invert = False
                    if mask_pattern == 0:
                        invert = ((r + c) % 2 == 0)
                    elif mask_pattern == 1:
                        invert = (r % 2 == 0)
                    elif mask_pattern == 2:
                        invert = (c % 3 == 0)
                    elif mask_pattern == 3:
                        invert = ((r + c) % 3 == 0)
                    elif mask_pattern == 4:
                        invert = ((r // 2 + c // 3) % 2 == 0)
                    elif mask_pattern == 5:
                        invert = (((r * c) % 2) + ((r * c) % 3) == 0)
                    elif mask_pattern == 6:
                        invert = ((((r * c) % 2) + ((r * c) % 3)) % 2 == 0)
                    elif mask_pattern == 7:
                        invert = ((((r + c) % 2) + ((r * c) % 3)) % 2 == 0)
                        
                    if invert:
                        val ^= 1
                    bitstream.append(val)
                    bit_idx += 1
        col -= 2
        upward = not upward
        
    # Map bit erasures to byte erasures
    erased_bytes = sorted(list(set(b // 8 for b in erasure_bits)))
    return bitstream, erased_bytes

def parse_byte_mode_message(data_bytes: list[int]) -> str:
    """Parse raw bytes starting with a Byte Mode header (0100) and character length indicator."""
    bits = []
    for b in data_bytes:
        for i in range(7, -1, -1):
            bits.append((b >> i) & 1)
            
    mode = 0
    for i in range(4):
        mode = (mode << 1) | bits[i]
        
    if mode != 4:
        raise ValueError(f"Unsupported mode indicator: {bin(mode)}. Only Byte Mode is supported.")
        
    char_len = 0
    for i in range(4, 12):
        char_len = (char_len << 1) | bits[i]
        
    msg_bytes = []
    bit_offset = 12
    for _ in range(char_len):
        if bit_offset + 8 > len(bits):
            break
        byte_val = 0
        for i in range(8):
            byte_val = (byte_val << 1) | bits[bit_offset + i]
        msg_bytes.append(byte_val)
        bit_offset += 8
        
    return bytes(msg_bytes).decode('utf-8', errors='replace')

def decode_qr_grid(grid: list[list[int]], types: list[list[str]], erase_locs: list[int] = None) -> tuple[str, list[int], list[int]]:
    """
    Decodes the given QR matrix:
    1. Read and decode format info.
    2. Extract and unmask data bits (collecting dynamic erasures from scratched modules).
    3. Convert to bytes.
    4. Apply Reed-Solomon correction.
    5. Parse mode & data to recover original text.
    """
    size = len(grid)
    version = 1 if size == 21 else 2
    
    left_bits, right_bits = extract_format_bits(grid)
    ec_level, mask_pattern = decode_format_info(left_bits, right_bits)
    
    total_data_ecc_bytes = 26 if version == 1 else 44
    data_bytes_cap, ecc_bytes, _ = QR_SPECS[(version, ec_level)]
    
    unmasked_bits, dynamic_erasures = extract_and_unmask_data_bits(grid, types, mask_pattern)
    
    combined_erasures = sorted(list(set((erase_locs or []) + dynamic_erasures)))
    # Only keep erasures that fit within the byte boundary limit
    combined_erasures = [e for e in combined_erasures if e < total_data_ecc_bytes]
    
    raw_bytes = bits_to_bytes(unmasked_bits)
    raw_bytes = raw_bytes[:total_data_ecc_bytes]
    
    corrected_bytes, err_pos, magnitudes = rs_correct_msg(raw_bytes, ecc_bytes, combined_erasures)
    
    data_bytes = corrected_bytes[:data_bytes_cap]
    decoded_text = parse_byte_mode_message(data_bytes)
    
    return decoded_text, err_pos, magnitudes

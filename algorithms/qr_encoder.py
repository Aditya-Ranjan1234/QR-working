# qr_encoder.py
# Pure Python QR Code Generator (Supports Version 1 & 2, Byte Mode, Levels L, M, Q, H)

from .reedsolomon import rs_encode_msg

# QR Code Specifications
# format: (version, ec_level) -> (data_bytes, ecc_bytes, matrix_size)
QR_SPECS = {
    (1, 'L'): (19, 7, 21),
    (1, 'M'): (16, 10, 21),
    (1, 'Q'): (13, 13, 21),
    (1, 'H'): (9, 17, 21),
    
    (2, 'L'): (34, 10, 25),
    (2, 'M'): (28, 16, 25),
    (2, 'Q'): (22, 22, 25),
    (2, 'H'): (16, 28, 25),
}

EC_INDICATORS = {
    'L': 0b01,
    'M': 0b00,
    'Q': 0b11,
    'H': 0b10,
}

# Format protection polynomial mask
FORMAT_MASK = 0x5412 # 101010000010010

def get_best_version(text: str, ec_level: str) -> int:
    """Find the smallest QR version (1 or 2) that can fit the text."""
    byte_len = len(text.encode('utf-8'))
    # Version 1 capacity (in bytes)
    v1_data_cap = QR_SPECS[(1, ec_level)][0]
    # Header: 4 bits mode + 8 bits length = 12 bits = 1.5 bytes.
    if byte_len <= v1_data_cap - 2:
        return 1
    v2_data_cap = QR_SPECS[(2, ec_level)][0]
    if byte_len <= v2_data_cap - 2:
        return 2
    raise ValueError(f"Text too long for simple Version 1 or 2 encoder under Level {ec_level}.")

def encode_data_to_bits(text: str, version: int, ec_level: str) -> list[int]:
    """Convert user data into a stream of bits, with headers, padding, and termination."""
    data_bytes_cap, _, _ = QR_SPECS[(version, ec_level)]
    data_bits_cap = data_bytes_cap * 8
    
    encoded_bytes = text.encode('utf-8')
    bitstream = []
    
    # 1. Mode Indicator: Byte mode is 0100 (4 bits)
    bitstream.extend([0, 1, 0, 0])
    
    # 2. Character Count Indicator: 8 bits for Version 1 and 2
    char_count = len(encoded_bytes)
    for i in range(7, -1, -1):
        bitstream.append((char_count >> i) & 1)
        
    # 3. Add character data
    for byte in encoded_bytes:
        for i in range(7, -1, -1):
            bitstream.append((byte >> i) & 1)
            
    # 4. Terminator: append up to 4 zero bits if capacity allows
    terminator_len = min(4, data_bits_cap - len(bitstream))
    bitstream.extend([0] * terminator_len)
    
    # 5. Round to nearest byte boundary
    while len(bitstream) % 8 != 0:
        bitstream.append(0)
        
    # 6. Stuff padding bytes (0xEC, 0x11 alternately) until capacity is reached
    padding_patterns = [0xEC, 0x11]
    pad_idx = 0
    while len(bitstream) < data_bits_cap:
        pad_byte = padding_patterns[pad_idx]
        pad_idx = 1 - pad_idx
        for i in range(7, -1, -1):
            bitstream.append((pad_byte >> i) & 1)
            
    return bitstream

def bits_to_bytes(bitstream: list[int]) -> list[int]:
    """Convert bitstream (list of 0s and 1s) to integer bytes."""
    bytes_list = []
    for i in range(0, len(bitstream), 8):
        byte_val = 0
        for bit in bitstream[i:i+8]:
            byte_val = (byte_val << 1) | bit
        bytes_list.append(byte_val)
    return bytes_list

def bytes_to_bits(bytes_list: list[int]) -> list[int]:
    """Convert list of bytes back to a bitstream."""
    bits = []
    for b in bytes_list:
        for i in range(7, -1, -1):
            bits.append((b >> i) & 1)
    return bits

# --- Format Information generator ---
def get_format_info(ec_level: str, mask_pattern: int) -> list[int]:
    """
    Calculate the 15-bit format info using generator polynomial G(x) = x^10 + x^8 + x^5 + x^4 + x^2 + x + 1 (0x537).
    """
    data = (EC_INDICATORS[ec_level] << 3) | mask_pattern
    # Calculate BCH error correction bits (10 bits)
    # data << 10
    rem = data << 10
    gen = 0x537 # generator polynomial binary representation
    for i in range(4, -1, -1):
        if rem & (1 << (10 + i)):
            rem ^= (gen << i)
            
    format_val = (data << 10) | rem
    # XOR with format mask
    masked_format = format_val ^ FORMAT_MASK
    
    # Return as 15 bits
    return [(masked_format >> i) & 1 for i in range(14, -1, -1)]

# --- QR Matrix Layout ---
class QRMatrix:
    def __init__(self, version: int):
        self.version = version
        self.size = QR_SPECS[(version, 'L')][2]
        # Grid values: None = empty, 0 = light, 1 = dark
        self.grid = [[None for _ in range(self.size)] for _ in range(self.size)]
        # Type grid: tracks what category each module belongs to (for coloring/visualization)
        self.types = [[ 'empty' for _ in range(self.size)] for _ in range(self.size)]
        
        self._place_static_patterns()
        
    def _place_static_patterns(self):
        # 1. Finder patterns at Top-Left, Top-Right, Bottom-Left
        self._place_finder(0, 0)
        self._place_finder(0, self.size - 7)
        self._place_finder(self.size - 7, 0)
        
        # 2. Separators around finders
        self._place_separators()
        
        # 3. Timing patterns (Row 6 and Column 6)
        for i in range(8, self.size - 8):
            # Alternating dark/light starting with dark (1) at index 8
            val = 1 if i % 2 == 0 else 0
            if self.grid[6][i] is None:
                self.grid[6][i] = val
                self.types[6][i] = 'timing'
            if self.grid[i][6] is None:
                self.grid[i][6] = val
                self.types[i][6] = 'timing'
                
        # 4. Alignment patterns
        if self.version == 2:
            self._place_alignment(18, 18)
            
        # 5. Dark module (always at column 8, row 4*version + 9)
        dark_row = 4 * self.version + 9
        self.grid[dark_row][8] = 1
        self.types[dark_row][8] = 'dark_module'
        
        # 6. Format info placeholders (around finders)
        # We reserve these modules so the data placement skips them
        self._reserve_format_info()
        
    def _place_finder(self, r: int, c: int):
        for dr in range(7):
            for dc in range(7):
                # Outer border (1), inner space (0), inner square (1)
                is_border = (dr == 0 or dr == 6 or dc == 0 or dc == 6)
                is_inner_square = (2 <= dr <= 4 and 2 <= dc <= 4)
                val = 1 if (is_border or is_inner_square) else 0
                self.grid[r + dr][c + dc] = val
                self.types[r + dr][c + dc] = 'finder'
                
    def _place_separators(self):
        # Top-Left separator
        for i in range(8):
            self._set_if_empty(7, i, 0, 'separator')
            self._set_if_empty(i, 7, 0, 'separator')
        # Top-Right separator
        for i in range(8):
            self._set_if_empty(7, self.size - 1 - i, 0, 'separator')
            self._set_if_empty(i, self.size - 8, 0, 'separator')
        # Bottom-Left separator
        for i in range(8):
            self._set_if_empty(self.size - 8, i, 0, 'separator')
            self._set_if_empty(self.size - 1 - i, 7, 0, 'separator')
            
    def _place_alignment(self, r: int, c: int):
        # 5x5 pattern centered at (r, c)
        for dr in range(-2, 3):
            for dc in range(-2, 3):
                # Outer border (1), inner space (0), center (1)
                is_border = (abs(dr) == 2 or abs(dc) == 2)
                is_center = (dr == 0 and dc == 0)
                val = 1 if (is_border or is_center) else 0
                self.grid[r + dr][c + dc] = val
                self.types[r + dr][c + dc] = 'alignment'
                
    def _reserve_format_info(self):
        # Left finder format area
        for i in range(9):
            if self.grid[i][8] is None:
                self.grid[i][8] = 0
                self.types[i][8] = 'format'
        for i in range(8):
            if self.grid[8][i] is None:
                self.grid[8][i] = 0
                self.types[8][i] = 'format'
                
        # Right/Bottom finder format area
        for i in range(8):
            r = self.size - 1 - i
            if self.grid[r][8] is None:
                self.grid[r][8] = 0
                self.types[r][8] = 'format'
        for i in range(8):
            c = self.size - 8 + i
            if self.grid[8][c] is None:
                self.grid[8][c] = 0
                self.types[8][c] = 'format'
                
    def _set_if_empty(self, r: int, c: int, val: int, t: str):
        if 0 <= r < self.size and 0 <= c < self.size:
            if self.grid[r][c] is None:
                self.grid[r][c] = val
                self.types[r][c] = t

    def place_data_bits(self, bitstream: list[int]):
        """Place data and ECC bitstream in empty modules using standard zig-zag path."""
        bit_idx = 0
        bit_len = len(bitstream)
        
        # Standard snake traversal: from bottom-right towards left, column-by-column in pairs
        col = self.size - 1
        upward = True
        
        while col > 0:
            # Skip vertical timing pattern column at index 6
            if col == 6:
                col -= 1
                
            rows = range(self.size - 1, -1, -1) if upward else range(self.size)
            for r in rows:
                # Place in two columns: col and col - 1
                for c in (col, col - 1):
                    if self.grid[r][c] is None or self.types[r][c] == 'empty':
                        bit_val = bitstream[bit_idx] if bit_idx < bit_len else 0
                        self.grid[r][c] = bit_val
                        self.types[r][c] = 'data'
                        bit_idx += 1
            col -= 2
            upward = not upward

    def apply_mask(self, mask_pattern: int):
        """Apply selected mask pattern to the data modules."""
        for r in range(self.size):
            for c in range(self.size):
                if self.types[r][c] == 'data':
                    # Determine if this coordinate gets inverted
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
                        self.grid[r][c] ^= 1

    def place_format_info(self, format_bits: list[int]):
        """Write the 15-bit format info into the reserved format modules."""
        # Left finder format area
        # format_bits: 15 bits, from most significant (index 0) to least (index 14)
        # Sequence on left/top: (0,8), (1,8), (2,8), (3,8), (4,8), (5,8), (7,8), (8,8), (8,7), (8,5), (8,4), (8,3), (8,2), (8,1), (8,0)
        # Let's map it safely:
        left_coords = [
            (0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (7, 8), (8, 8),
            (8, 7), (8, 5), (8, 4), (8, 3), (8, 2), (8, 1), (8, 0)
        ]
        
        # Right/bottom format area:
        # Sequence: (8, size-1) to (8, size-8) for bits 0..7
        # (size-7, 8) to (size-1, 8) for bits 8..14
        right_coords = []
        for i in range(8):
            right_coords.append((8, self.size - 1 - i))
        for i in range(7):
            right_coords.append((self.size - 7 + i, 8))
            
        for i in range(15):
            r_l, c_l = left_coords[i]
            self.grid[r_l][c_l] = format_bits[i]
            self.types[r_l][c_l] = 'format'
            
            r_r, c_r = right_coords[i]
            self.grid[r_r][c_r] = format_bits[i]
            self.types[r_r][c_r] = 'format'

    # --- Penalty Calculations (Rule 1-4) ---
    def calculate_penalty(self) -> int:
        penalty = 0
        
        # Rule 1: 5 or more consecutive modules of the same color in row/col
        # Row-wise
        for r in range(self.size):
            run = 1
            for c in range(1, self.size):
                if self.grid[r][c] == self.grid[r][c-1]:
                    run += 1
                else:
                    if run >= 5:
                        penalty += 3 + (run - 5)
                    run = 1
            if run >= 5:
                penalty += 3 + (run - 5)
                
        # Col-wise
        for c in range(self.size):
            run = 1
            for r in range(1, self.size):
                if self.grid[r][c] == self.grid[r-1][c]:
                    run += 1
                else:
                    if run >= 5:
                        penalty += 3 + (run - 5)
                    run = 1
            if run >= 5:
                penalty += 3 + (run - 5)
                
        # Rule 2: 2x2 blocks of the same color
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                v = self.grid[r][c]
                if v == self.grid[r+1][c] == self.grid[r][c+1] == self.grid[r+1][c+1]:
                    penalty += 3
                    
        # Rule 3: Dark-Light-Dark-Dark-Dark-Light-Dark-Light-Light-Light-Light (1-0-1-1-1-0-1-0-0-0-0)
        # or reversed (0-0-0-0-1-0-1-1-1-0-1) in row or column
        p1 = [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0]
        p2 = [0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1]
        
        for r in range(self.size):
            for c in range(self.size - 10):
                slice_row = self.grid[r][c:c+11]
                if slice_row == p1 or slice_row == p2:
                    penalty += 40
                    
        for c in range(self.size):
            for r in range(self.size - 10):
                slice_col = [self.grid[r+i][c] for i in range(11)]
                if slice_col == p1 or slice_col == p2:
                    penalty += 40
                    
        # Rule 4: Proportion of dark modules
        total_modules = self.size * self.size
        dark_modules = sum(sum(row) for row in self.grid)
        percentage = (dark_modules / total_modules) * 100
        deviation = abs(percentage - 50)
        penalty += int(deviation // 5) * 10
        
        return penalty

    def calculate_penalty_detailed(self) -> dict:
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        
        # Rule 1
        for r in range(self.size):
            run = 1
            for c in range(1, self.size):
                if self.grid[r][c] == self.grid[r][c-1]:
                    run += 1
                else:
                    if run >= 5:
                        n1 += 3 + (run - 5)
                    run = 1
            if run >= 5:
                n1 += 3 + (run - 5)
                
        for c in range(self.size):
            run = 1
            for r in range(1, self.size):
                if self.grid[r][c] == self.grid[r-1][c]:
                    run += 1
                else:
                    if run >= 5:
                        n1 += 3 + (run - 5)
                    run = 1
            if run >= 5:
                n1 += 3 + (run - 5)
                
        # Rule 2
        for r in range(self.size - 1):
            for c in range(self.size - 1):
                v = self.grid[r][c]
                if v == self.grid[r+1][c] == self.grid[r][c+1] == self.grid[r+1][c+1]:
                    n2 += 3
                    
        # Rule 3
        p1 = [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0]
        p2 = [0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1]
        
        for r in range(self.size):
            for c in range(self.size - 10):
                slice_row = self.grid[r][c:c+11]
                if slice_row == p1 or slice_row == p2:
                    n3 += 40
                    
        for c in range(self.size):
            for r in range(self.size - 10):
                slice_col = [self.grid[r+i][c] for i in range(11)]
                if slice_col == p1 or slice_col == p2:
                    n3 += 40
                    
        # Rule 4
        total_modules = self.size * self.size
        dark_modules = sum(sum(row) for row in self.grid)
        percentage = (dark_modules / total_modules) * 100
        deviation = abs(percentage - 50)
        n4 += int(deviation // 5) * 10
        
        return {
            "n1": n1,
            "n2": n2,
            "n3": n3,
            "n4": n4,
            "total": n1 + n2 + n3 + n4
        }

def generate_qr_data_and_ecc(text: str, version: int, ec_level: str) -> list[int]:
    """Full data & ECC bytes pipeline."""
    # 1. Get raw bitstream
    data_bits = encode_data_to_bits(text, version, ec_level)
    # 2. Group into data bytes
    data_bytes = bits_to_bytes(data_bits)
    # 3. Call Reed-Solomon encoder to get ECC
    _, ecc_bytes, _ = QR_SPECS[(version, ec_level)]
    full_bytes = rs_encode_msg(data_bytes, ecc_bytes)
    # 4. Turn full bytes (data + ECC) back to bitstream for matrix placement
    return bytes_to_bits(full_bytes)

def create_qr_code(text: str, ec_level: str = 'M') -> tuple[QRMatrix, int, list[int]]:
    """Determine version, generate, mask, and return the final matrix, selected mask, and byte list."""
    version = get_best_version(text, ec_level)
    _, ecc_bytes, _ = QR_SPECS[(version, ec_level)]
    
    # Generate data + ecc bitstream
    full_bitstream = generate_qr_data_and_ecc(text, version, ec_level)
    full_bytes = bits_to_bytes(full_bitstream)
    
    # Try all 8 masks and pick the best one
    best_mask = 0
    best_penalty = 99999999
    best_matrix = None
    
    for mask in range(8):
        mat = QRMatrix(version)
        mat.place_data_bits(full_bitstream)
        mat.apply_mask(mask)
        format_bits = get_format_info(ec_level, mask)
        mat.place_format_info(format_bits)
        
        penalty = mat.calculate_penalty()
        if penalty < best_penalty:
            best_penalty = penalty
            best_mask = mask
            best_matrix = mat
            
    return best_matrix, best_mask, full_bytes

def create_qr_code_verbose(text: str, ec_level: str = 'M') -> tuple[QRMatrix, int, list[int], dict]:
    """Verbose creator returning final QR details along with all 8 mask penalty score details."""
    version = get_best_version(text, ec_level)
    _, ecc_bytes, _ = QR_SPECS[(version, ec_level)]
    full_bitstream = generate_qr_data_and_ecc(text, version, ec_level)
    full_bytes = bits_to_bytes(full_bitstream)
    
    best_mask = 0
    best_penalty = 99999999
    best_matrix = None
    mask_details = {}
    
    for mask in range(8):
        mat = QRMatrix(version)
        mat.place_data_bits(full_bitstream)
        mat.apply_mask(mask)
        format_bits = get_format_info(ec_level, mask)
        mat.place_format_info(format_bits)
        
        details = mat.calculate_penalty_detailed()
        mask_details[mask] = details
        
        if details["total"] < best_penalty:
            best_penalty = details["total"]
            best_mask = mask
            best_matrix = mat
            
    return best_matrix, best_mask, full_bytes, mask_details

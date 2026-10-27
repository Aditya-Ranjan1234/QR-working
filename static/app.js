// app.js
// Complete Serverless Client-Side Galois Field, Reed-Solomon and QR Code Engine

// ==========================================
// 1. GALOIS FIELD GF(256) ARITHMETIC
// ==========================================
const PRIMITIVE_POLYNOMIAL = 285; // x^8 + x^4 + x^3 + x^2 + 1
let gf_exp = new Array(512).fill(0);
let gf_log = new Array(256).fill(0);

function init_gf256_tables() {
    let x = 1;
    for (let i = 0; i < 255; i++) {
        gf_exp[i] = x;
        gf_log[x] = i;
        x <<= 1;
        if (x & 0x100) {
            x ^= PRIMITIVE_POLYNOMIAL;
        }
    }
    for (let i = 255; i < 512; i++) {
        gf_exp[i] = gf_exp[i - 255];
    }
}
init_gf256_tables();

function gf_add(a, b) { return a ^ b; }
function gf_sub(a, b) { return a ^ b; }

function gf_mul(a, b) {
    if (a === 0 || b === 0) return 0;
    return gf_exp[gf_log[a] + gf_log[b]];
}

function gf_inv(a) {
    if (a === 0) throw new Error("GF(256) division by zero.");
    return gf_exp[255 - gf_log[a]];
}

function gf_div(a, b) {
    if (b === 0) throw new Error("GF(256) division by zero.");
    if (a === 0) return 0;
    return gf_exp[gf_log[a] + 255 - gf_log[b]];
}

function gf_pow(a, power) {
    if (a === 0) return power > 0 ? 0 : 1;
    return gf_exp[(gf_log[a] * power) % 255];
}

// Format byte as polynomial string helper
function byte_to_poly_str(val) {
    if (val === 0) return "0";
    let terms = [];
    for (let i = 7; i >= 0; i--) {
        if ((val >> i) & 1) {
            if (i === 0) terms.push("1");
            else if (i === 1) terms.push("x");
            else terms.push(`x^${i}`);
        }
    }
    return terms.join(" + ");
}

// ==========================================
// 2. POLYNOMIAL ARITHMETIC OVER GF(256)
// ==========================================
function poly_add(p1, p2) {
    let len = Math.max(p1.length, p2.length);
    let r = new Array(len).fill(0);
    for (let i = 0; i < p1.length; i++) {
        r[i + len - p1.length] ^= p1[i];
    }
    for (let i = 0; i < p2.length; i++) {
        r[i + len - p2.length] ^= p2[i];
    }
    return r;
}

function poly_mul(p1, p2) {
    let r = new Array(p1.length + p2.length - 1).fill(0);
    for (let i = 0; i < p1.length; i++) {
        for (let j = 0; j < p2.length; j++) {
            r[i + j] ^= gf_mul(p1[i], p2[j]);
        }
    }
    return r;
}

function poly_eval(p, x) {
    let y = p[0];
    for (let i = 1; i < p.length; i++) {
        y = gf_add(gf_mul(y, x), p[i]);
    }
    return y;
}

function poly_div(num, den) {
    let r = [...num];
    for (let i = 0; i < num.length - den.length + 1; i++) {
        let coef = r[i];
        if (coef !== 0) {
            for (let j = 1; j < den.length; j++) {
                r[i + j] ^= gf_mul(den[j], coef);
            }
        }
    }
    let separator = num.length - den.length + 1;
    return [r.slice(0, separator), r.slice(separator)];
}

// Polynomial helpers with lowest-degree first ordering for Reed-Solomon
function gf_poly_mul(p, q) {
    let r = new Array(p.length + q.length - 1).fill(0);
    for (let j = 0; j < q.length; j++) {
        for (let i = 0; i < p.length; i++) {
            r[i + j] ^= gf_mul(p[i], q[j]);
        }
    }
    return r;
}

function gf_poly_eval(p, x) {
    let y = 0;
    let xp = 1;
    for (let coef of p) {
        y ^= gf_mul(coef, xp);
        xp = gf_mul(xp, x);
    }
    return y;
}

// ==========================================
// 3. REED-SOLOMON CODEC (Client-Side)
// ==========================================
function rs_generator_poly(nsym) {
    let g = [1];
    for (let i = 0; i < nsym; i++) {
        g = gf_poly_mul(g, [gf_exp[i], 1]);
    }
    return g;
}

function rs_encode_msg(msg, nsym) {
    let gen = rs_generator_poly(nsym);
    let msg_out = new Array(nsym).fill(0).concat([...msg].reverse());
    
    for (let i = msg.length - 1; i >= 0; i--) {
        let coef = msg_out[i + nsym];
        if (coef !== 0) {
            for (let j = 0; j < gen.length - 1; j++) {
                msg_out[i + j] ^= gf_mul(gen[j], coef);
            }
        }
    }
    let parity = msg_out.slice(0, nsym);
    return [...msg].concat(parity.reverse());
}

function rs_calc_syndromes(msg, nsym) {
    let msg_rev = [...msg].reverse();
    let synd = new Array(nsym).fill(0);
    for (let i = 0; i < nsym; i++) {
        synd[i] = gf_poly_eval(msg_rev, gf_exp[i]);
    }
    return [0].concat(synd);
}

function rs_find_error_locator_verbose(synd, nsym, erase_loc = []) {
    let lam = [1];
    for (let el of erase_loc) {
        let X = gf_exp[el];
        lam = gf_poly_mul(lam, [X, 1]);
    }
    
    let B = [...lam];
    let L = erase_loc.length;
    let steps = [];
    
    steps.push({
        r: 0, d: 0, L: L,
        lam: [...lam], B: [...B],
        description: erase_loc.length > 0 ? "Initial state (erasures locator polynomial initialized)" : "Initial state"
    });
    
    for (let r = L + 1; r <= nsym; r++) {
        let d = synd[r];
        for (let j = 1; j < lam.length; j++) {
            if (r - j > 0) {
                d ^= gf_mul(lam[j], synd[r - j]);
            }
        }
        
        B = [0].concat(B);
        let desc = `r=${r}: calculated discrepancy d=${d}. `;
        
        if (d !== 0) {
            if (2 * L < r + erase_loc.length) {
                let old_lam = [...lam];
                let scaled_B = B.map(x => gf_mul(x, d));
                if (lam.length < scaled_B.length) {
                    lam = lam.concat(new Array(scaled_B.length - lam.length).fill(0));
                }
                for (let j = 0; j < scaled_B.length; j++) {
                    lam[j] ^= scaled_B[j];
                }
                L = r + erase_loc.length - L;
                let inv_d = gf_inv(d);
                B = old_lam.map(x => gf_mul(x, inv_d));
                desc += `d != 0 and 2*L < r + erasures. Lambda updated. L updated to ${L}. B updated to old_Lambda / d.`;
            } else {
                let scaled_B = B.map(x => gf_mul(x, d));
                if (lam.length < scaled_B.length) {
                    lam = lam.concat(new Array(scaled_B.length - lam.length).fill(0));
                }
                for (let j = 0; j < scaled_B.length; j++) {
                    lam[j] ^= scaled_B[j];
                }
                desc += "d != 0 but 2*L >= r + erasures. Lambda updated. L unchanged. B shifted.";
            }
        } else {
            desc += "d == 0. No update to Lambda. L unchanged. B shifted.";
        }
        
        let lam_clean = [...lam];
        while (lam_clean.length > 0 && lam_clean[lam_clean.length - 1] === 0) {
            lam_clean.pop();
        }
        
        steps.push({
            r: r, d: d, L: L,
            lam: lam_clean, B: [...B],
            description: desc
        });
    }
    
    while (lam.length > 0 && lam[lam.length - 1] === 0) {
        lam.pop();
    }
    
    return [lam, steps];
}

function rs_find_errors(err_loc, msg_len) {
    let errors = [];
    for (let i = 0; i < msg_len; i++) {
        let Xi_inv = gf_exp[255 - i];
        let val = gf_poly_eval(err_loc, Xi_inv);
        if (val === 0) {
            errors.push(msg_len - 1 - i);
        }
    }
    return errors;
}

function rs_forney_magnitudes(synd, err_loc, err_pos, msg_len) {
    let nsym = synd.length - 1;
    let S = synd.slice(1);
    let omega = gf_poly_mul(S, err_loc).slice(0, nsym);
    
    let lam_prime = [];
    for (let i = 1; i < err_loc.length; i++) {
        lam_prime.push(i % 2 === 1 ? err_loc[i] : 0);
    }
    
    let magnitudes = [];
    for (let pos of err_pos) {
        let Xi = gf_exp[msg_len - 1 - pos];
        let Xi_inv = gf_inv(Xi);
        
        let num = gf_poly_eval(omega, Xi_inv);
        let denom = gf_poly_eval(lam_prime, Xi_inv);
        
        let mag = gf_mul(Xi, gf_div(num, denom));
        magnitudes.push(mag);
    }
    return magnitudes;
}

function rs_correct_msg(msg, nsym, erase_loc = []) {
    let corrected = [...msg];
    let synd = rs_calc_syndromes(msg, nsym);
    if (Math.max(...synd) === 0) {
        return [corrected, [], []];
    }
    
    let erase_loc_end = erase_loc.map(e => msg.length - 1 - e);
    let [err_loc, ] = rs_find_error_locator_verbose(synd, nsym, erase_loc_end);
    let err_pos = rs_find_errors(err_loc, msg.length);
    
    if (err_pos.length === 0) {
        throw new Error("Could not locate errors (too many errors).");
    }
    if (err_pos.length * 2 - erase_loc.length > nsym) {
        throw new Error("Too many errors/erasures to correct.");
    }
    
    let magnitudes = rs_forney_magnitudes(synd, err_loc, err_pos, msg.length);
    for (let i = 0; i < err_pos.length; i++) {
        corrected[err_pos[i]] ^= magnitudes[i];
    }
    
    let final_synd = rs_calc_syndromes(corrected, nsym);
    if (Math.max(...final_synd) !== 0) {
        throw new Error("Decoder failed to repair all errors.");
    }
    
    return [corrected, err_pos, magnitudes];
}

// ==========================================
// 4. QR CODE ENCODER ENGINE
// ==========================================
const QR_SPECS = {
    "1_L": [19, 7, 21], "1_M": [16, 10, 21], "1_Q": [13, 13, 21], "1_H": [9, 17, 21],
    "2_L": [34, 10, 25], "2_M": [28, 16, 25], "2_Q": [22, 22, 25], "2_H": [16, 28, 25]
};

const EC_INDICATORS = { L: 1, M: 0, Q: 3, H: 2 };
const FORMAT_MASK = 0x5412;

function get_best_version(text, ec_level) {
    let byte_len = new TextEncoder().encode(text).length;
    let v1_cap = QR_SPECS[`1_${ec_level}`][0];
    if (byte_len <= v1_cap - 2) return 1;
    let v2_cap = QR_SPECS[`2_${ec_level}`][0];
    if (byte_len <= v2_cap - 2) return 2;
    throw new Error(`Text too long for Version 1 or 2 under Level ${ec_level}.`);
}

function encode_data_to_bits(text, version, ec_level) {
    let [data_bytes_cap, , ] = QR_SPECS[`${version}_${ec_level}`];
    let data_bits_cap = data_bytes_cap * 8;
    
    let bytes = new TextEncoder().encode(text);
    let bits = [];
    
    // Mode Indicator: Byte mode is 0100
    bits.push(0, 1, 0, 0);
    
    // Character count indicator (8 bits)
    let char_count = bytes.length;
    for (let i = 7; i >= 0; i--) bits.push((char_count >> i) & 1);
    
    // Character data
    for (let b of bytes) {
        for (let i = 7; i >= 0; i--) bits.push((b >> i) & 1);
    }
    
    // Terminator
    let term_len = Math.min(4, data_bits_cap - bits.length);
    for (let i = 0; i < term_len; i++) bits.push(0);
    
    // Pad to 8-bit boundary
    while (bits.length % 8 !== 0) bits.push(0);
    
    // Stuff padding bytes 0xEC, 0x11
    let pads = [0xEC, 0x11];
    let pad_idx = 0;
    while (bits.length < data_bits_cap) {
        let pb = pads[pad_idx];
        pad_idx = 1 - pad_idx;
        for (let i = 7; i >= 0; i--) bits.push((pb >> i) & 1);
    }
    return bits;
}

function bits_to_bytes(bits) {
    let bytes = [];
    for (let i = 0; i < bits.length; i += 8) {
        let val = 0;
        for (let j = 0; j < 8; j++) val = (val << 1) | (bits[i + j] || 0);
        bytes.push(val);
    }
    return bytes;
}

function bytes_to_bits(bytes) {
    let bits = [];
    for (let b of bytes) {
        for (let i = 7; i >= 0; i--) bits.push((b >> i) & 1);
    }
    return bits;
}

function get_format_info(ec_level, mask) {
    let data = (EC_INDICATORS[ec_level] << 3) | mask;
    let rem = data << 10;
    let gen = 0x537;
    for (let i = 4; i >= 0; i--) {
        if (rem & (1 << (10 + i))) {
            rem ^= (gen << i);
        }
    }
    let format_val = ((data << 10) | rem) ^ FORMAT_MASK;
    let bits = [];
    for (let i = 14; i >= 0; i--) bits.push((format_val >> i) & 1);
    return bits;
}

class QRMatrix {
    constructor(version) {
        this.version = version;
        this.size = QR_SPECS[`${version}_L`][2];
        this.grid = Array.from({ length: this.size }, () => new Array(this.size).fill(null));
        this.types = Array.from({ length: this.size }, () => new Array(this.size).fill('empty'));
        this.place_static_patterns();
    }
    
    place_static_patterns() {
        this.place_finder(0, 0);
        this.place_finder(0, this.size - 7);
        this.place_finder(this.size - 7, 0);
        this.place_separators();
        
        for (let i = 8; i < this.size - 8; i++) {
            let val = i % 2 === 0 ? 1 : 0;
            if (this.grid[6][i] === null) {
                this.grid[6][i] = val;
                this.types[6][i] = 'timing';
            }
            if (this.grid[i][6] === null) {
                this.grid[i][6] = val;
                this.types[i][6] = 'timing';
            }
        }
        if (this.version === 2) {
            this.place_alignment(18, 18);
        }
        let dark_row = 4 * this.version + 9;
        this.grid[dark_row][8] = 1;
        this.types[dark_row][8] = 'dark_module';
        this.reserve_format_info();
    }
    
    place_finder(r, c) {
        for (let dr = 0; dr < 7; dr++) {
            for (let dc = 0; dc < 7; dc++) {
                let border = (dr === 0 || dr === 6 || dc === 0 || dc === 6);
                let center = (dr >= 2 && dr <= 4 && dc >= 2 && dc <= 4);
                this.grid[r + dr][c + dc] = (border || center) ? 1 : 0;
                this.types[r + dr][c + dc] = 'finder';
            }
        }
    }
    
    place_separators() {
        let setSep = (r, c) => {
            if (r >= 0 && r < this.size && c >= 0 && c < this.size) {
                if (this.grid[r][c] === null) {
                    this.grid[r][c] = 0;
                    this.types[r][c] = 'separator';
                }
            }
        };
        for (let i = 0; i < 8; i++) {
            setSep(7, i); setSep(i, 7);
            setSep(7, this.size - 1 - i); setSep(i, this.size - 8);
            setSep(this.size - 8, i); setSep(this.size - 1 - i, 7);
        }
    }
    
    place_alignment(r, c) {
        for (let dr = -2; dr <= 2; dr++) {
            for (let dc = -2; dc <= 2; dc++) {
                let border = (Math.abs(dr) === 2 || Math.abs(dc) === 2);
                let center = (dr === 0 && dc === 0);
                this.grid[r + dr][c + dc] = (border || center) ? 1 : 0;
                this.types[r + dr][c + dc] = 'alignment';
            }
        }
    }
    
    reserve_format_info() {
        for (let i = 0; i < 9; i++) {
            if (this.grid[i][8] === null) { this.grid[i][8] = 0; this.types[i][8] = 'format'; }
        }
        for (let i = 0; i < 8; i++) {
            if (this.grid[8][i] === null) { this.grid[8][i] = 0; this.types[8][i] = 'format'; }
        }
        for (let i = 0; i < 8; i++) {
            let r = this.size - 1 - i;
            if (this.grid[r][8] === null) { this.grid[r][8] = 0; this.types[r][8] = 'format'; }
        }
        for (let i = 0; i < 8; i++) {
            let c = this.size - 8 + i;
            if (this.grid[8][c] === null) { this.grid[8][c] = 0; this.types[8][c] = 'format'; }
        }
    }
    
    place_data_bits(bits) {
        let bit_idx = 0;
        let col = this.size - 1;
        let upward = true;
        while (col > 0) {
            if (col === 6) col--;
            let rows = upward ? Array.from({length: this.size}, (_, i) => this.size - 1 - i) : Array.from({length: this.size}, (_, i) => i);
            for (let r of rows) {
                for (let c of [col, col - 1]) {
                    if (this.grid[r][c] === null || this.types[r][c] === 'empty') {
                        this.grid[r][c] = bit_idx < bits.length ? bits[bit_idx] : 0;
                        this.types[r][c] = 'data';
                        bit_idx++;
                    }
                }
            }
            col -= 2;
            upward = !upward;
        }
    }
    
    apply_mask(mask) {
        for (let r = 0; r < this.size; r++) {
            for (let c = 0; c < this.size; c++) {
                if (this.types[r][c] === 'data') {
                    let invert = false;
                    if (mask === 0) invert = ((r + c) % 2 === 0);
                    else if (mask === 1) invert = (r % 2 === 0);
                    else if (mask === 2) invert = (c % 3 === 0);
                    else if (mask === 3) invert = ((r + c) % 3 === 0);
                    else if (mask === 4) invert = ((Math.floor(r / 2) + Math.floor(c / 3)) % 2 === 0);
                    else if (mask === 5) invert = (((r * c) % 2) + ((r * c) % 3) === 0);
                    else if (mask === 6) invert = ((((r * c) % 2) + ((r * c) % 3)) % 2 === 0);
                    else if (mask === 7) invert = ((((r + c) % 2) + ((r * c) % 3)) % 2 === 0);
                    
                    if (invert) this.grid[r][c] ^= 1;
                }
            }
        }
    }
    
    place_format_info(bits) {
        let left_coords = [
            [0, 8], [1, 8], [2, 8], [3, 8], [4, 8], [5, 8], [7, 8], [8, 8],
            [8, 7], [8, 5], [8, 4], [8, 3], [8, 2], [8, 1], [8, 0]
        ];
        let right_coords = [];
        for (let i = 0; i < 8; i++) right_coords.push([8, this.size - 1 - i]);
        for (let i = 0; i < 7; i++) right_coords.push([this.size - 7 + i, 8]);
        
        for (let i = 0; i < 15; i++) {
            let [rl, cl] = left_coords[i];
            this.grid[rl][cl] = bits[i];
            this.types[rl][cl] = 'format';
            
            let [rr, cr] = right_coords[i];
            this.grid[rr][cr] = bits[i];
            this.types[rr][cr] = 'format';
        }
    }
    
    calculate_penalty_detailed() {
        let n1 = 0, n2 = 0, n3 = 0, n4 = 0;
        // Rule 1: runs of 5 or more
        for (let r = 0; r < this.size; r++) {
            let run = 1;
            for (let c = 1; c < this.size; c++) {
                if (this.grid[r][c] === this.grid[r][c - 1]) run++;
                else { if (run >= 5) n1 += 3 + (run - 5); run = 1; }
            }
            if (run >= 5) n1 += 3 + (run - 5);
        }
        for (let c = 0; c < this.size; c++) {
            let run = 1;
            for (let r = 1; r < this.size; r++) {
                if (this.grid[r][c] === this.grid[r - 1][c]) run++;
                else { if (run >= 5) n1 += 3 + (run - 5); run = 1; }
            }
            if (run >= 5) n1 += 3 + (run - 5);
        }
        
        // Rule 2: 2x2 blocks
        for (let r = 0; r < this.size - 1; r++) {
            for (let c = 0; c < this.size - 1; c++) {
                let v = this.grid[r][c];
                if (v === this.grid[r + 1][c] && v === this.grid[r][c + 1] && v === this.grid[r + 1][c + 1]) {
                    n2 += 3;
                }
            }
        }
        
        // Rule 3: Finder pattern lookalikes
        let p1 = [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0];
        let p2 = [0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1];
        let matchArr = (a1, a2) => a1.length === a2.length && a1.every((v, idx) => v === a2[idx]);
        
        for (let r = 0; r < this.size; r++) {
            for (let c = 0; c < this.size - 10; c++) {
                let slice = this.grid[r].slice(c, c + 11);
                if (matchArr(slice, p1) || matchArr(slice, p2)) n3 += 40;
            }
        }
        for (let c = 0; c < this.size; c++) {
            for (let r = 0; r < this.size - 10; r++) {
                let slice = [];
                for (let i = 0; i < 11; i++) slice.push(this.grid[r + i][c]);
                if (matchArr(slice, p1) || matchArr(slice, p2)) n3 += 40;
            }
        }
        
        // Rule 4: Balance of dark/light modules
        let total = this.size * this.size;
        let dark = 0;
        for (let r = 0; r < this.size; r++) {
            for (let c = 0; c < this.size; c++) if (this.grid[r][c] === 1) dark++;
        }
        let percentage = (dark / total) * 100;
        let deviation = Math.abs(percentage - 50);
        n4 = Math.floor(deviation / 5) * 10;
        
        return { n1, n2, n3, n4, total: n1 + n2 + n3 + n4 };
    }
}

function generate_qr_data_and_ecc(text, version, ec_level) {
    let bits = encode_data_to_bits(text, version, ec_level);
    let bytes = bits_to_bytes(bits);
    let [, ecc_bytes_count, ] = QR_SPECS[`${version}_${ec_level}`];
    let full = rs_encode_msg(bytes, ecc_bytes_count);
    return bytes_to_bits(full);
}

function create_qr_code_verbose(text, ec_level) {
    let version = get_best_version(text, ec_level);
    let full_bits = generate_qr_data_and_ecc(text, version, ec_level);
    let full_bytes = bits_to_bytes(full_bits);
    
    let best_mask = 0;
    let best_penalty = Infinity;
    let best_matrix = null;
    let mask_details = {};
    
    for (let mask = 0; mask < 8; mask++) {
        let mat = new QRMatrix(version);
        mat.place_data_bits(full_bits);
        mat.apply_mask(mask);
        let format_bits = get_format_info(ec_level, mask);
        mat.place_format_info(format_bits);
        
        let details = mat.calculate_penalty_detailed();
        mask_details[mask] = details;
        
        if (details.total < best_penalty) {
            best_penalty = details.total;
            best_mask = mask;
            best_matrix = mat;
        }
    }
    return [best_matrix, best_mask, full_bytes, mask_details];
}

// ==========================================
// 5. QR CODE DECODER ENGINE
// ==========================================
function decode_format_info(left_bits, right_bits) {
    let best_ec = 'M';
    let best_mask = 0;
    let best_dist = Infinity;
    
    let matchArr = (a1, a2) => a1.every((v, i) => v === a2[i]);
    
    for (let ec of ['L', 'M', 'Q', 'H']) {
        for (let mask = 0; mask < 8; mask++) {
            let valid = get_format_info(ec, mask);
            
            let dist_l = left_bits.reduce((acc, v, i) => acc + (v !== valid[i] ? 1 : 0), 0);
            let dist_r = right_bits.reduce((acc, v, i) => acc + (v !== valid[i] ? 1 : 0), 0);
            
            let min_d = Math.min(dist_l, dist_r);
            if (min_d < best_dist) {
                best_dist = min_d;
                best_ec = ec;
                best_mask = mask;
            }
        }
    }
    return [best_ec, best_mask];
}

function extract_format_bits(grid) {
    let size = grid.length;
    let left_coords = [
        [0, 8], [1, 8], [2, 8], [3, 8], [4, 8], [5, 8], [7, 8], [8, 8],
        [8, 7], [8, 5], [8, 4], [8, 3], [8, 2], [8, 1], [8, 0]
    ];
    let right_coords = [];
    for (let i = 0; i < 8; i++) right_coords.push([8, size - 1 - i]);
    for (let i = 0; i < 7; i++) right_coords.push([size - 7 + i, 8]);
    
    let left = left_coords.map(([r, c]) => grid[r][c] === null ? 0 : grid[r][c]);
    let right = right_coords.map(([r, c]) => grid[r][c] === null ? 0 : grid[r][c]);
    
    return [left, right];
}

function extract_and_unmask_data_bits(grid, types, mask) {
    let size = grid.length;
    let bits = [];
    let erasure_bits = [];
    
    let col = size - 1;
    let upward = true;
    let bit_idx = 0;
    
    while (col > 0) {
        if (col === 6) col--;
        let rows = upward ? Array.from({length: size}, (_, i) => size - 1 - i) : Array.from({length: size}, (_, i) => i);
        for (let r of rows) {
            for (let c of [col, col - 1]) {
                if (types[r][c] === 'data') {
                    let val = grid[r][c];
                    if (val === -1) {
                        erasure_bits.push(bit_idx);
                        val = 0;
                    }
                    
                    let invert = false;
                    if (mask === 0) invert = ((r + c) % 2 === 0);
                    else if (mask === 1) invert = (r % 2 === 0);
                    else if (mask === 2) invert = (c % 3 === 0);
                    else if (mask === 3) invert = ((r + c) % 3 === 0);
                    else if (mask === 4) invert = ((Math.floor(r / 2) + Math.floor(c / 3)) % 2 === 0);
                    else if (mask === 5) invert = (((r * c) % 2) + ((r * c) % 3) === 0);
                    else if (mask === 6) invert = ((((r * c) % 2) + ((r * c) % 3)) % 2 === 0);
                    else if (mask === 7) invert = ((((r + c) % 2) + ((r * c) % 3)) % 2 === 0);
                    
                    if (invert) val ^= 1;
                    bits.push(val);
                    bit_idx++;
                }
            }
        }
        col -= 2;
        upward = !upward;
    }
    
    let erased_bytes = sortedUnique(erasure_bits.map(b => Math.floor(b / 8)));
    return [bits, erased_bytes];
}

function sortedUnique(arr) {
    return Array.from(new Set(arr)).sort((a, b) => a - b);
}

function parse_byte_mode_message(bytes) {
    let bits = [];
    for (let b of bytes) {
        for (let i = 7; i >= 0; i--) bits.push((b >> i) & 1);
    }
    
    let mode = 0;
    for (let i = 0; i < 4; i++) mode = (mode << 1) | bits[i];
    if (mode !== 4) throw new Error(`Unsupported QR mode indicator: ${mode}. Only Byte Mode (0100) is supported.`);
    
    let char_len = 0;
    for (let i = 4; i < 12; i++) char_len = (char_len << 1) | bits[i];
    
    let msg_bytes = [];
    let bit_offset = 12;
    for (let i = 0; i < char_len; i++) {
        if (bit_offset + 8 > bits.length) break;
        let val = 0;
        for (let j = 0; j < 8; j++) val = (val << 1) | bits[bit_offset + j];
        msg_bytes.push(val);
        bit_offset += 8;
    }
    return new TextDecoder().decode(new Uint8Array(msg_bytes));
}

function decode_qr_grid_verbose(grid, types, erased_idx = []) {
    let size = grid.length;
    let version = size === 21 ? 1 : 2;
    
    let [left, right] = extract_format_bits(grid);
    let [ec_level, mask] = decode_format_info(left, right);
    
    let total_bytes = version === 1 ? 26 : 44;
    let [data_cap, ecc_bytes_count, ] = QR_SPECS[`${version}_${ec_level}`];
    
    let [unmasked_bits, dynamic_erasures] = extract_and_unmask_data_bits(grid, types, mask);
    let raw_bytes = bits_to_bytes(unmasked_bits).slice(0, total_bytes);
    
    let combined_erasures = sortedUnique(erased_idx.concat(dynamic_erasures)).filter(e => e < total_bytes);
    let combined_erasures_end = combined_erasures.map(e => raw_bytes.length - 1 - e);
    
    let synd = rs_calc_syndromes(raw_bytes, ecc_bytes_count);
    let [err_loc, bm_steps] = rs_find_error_locator_verbose(synd, ecc_bytes_count, combined_erasures_end);
    
    let [corrected, err_pos, magnitudes] = rs_correct_msg(raw_bytes, ecc_bytes_count, combined_erasures);
    let parsed_text = parse_byte_mode_message(corrected.slice(0, data_cap));
    
    let logs = [
        `Version: ${version} | EC Level: ${ec_level} | Mask: ${mask}`,
        `Extracted ${raw_bytes.length} bytes from grid zig-zag layout.`,
        combined_erasures.length > 0 ? `Erasures at codeword indices: [${combined_erasures.join(', ')}]` : "No known erasures.",
        err_pos.length > 0 ? `Repaired ${err_pos.length} errors at indices: [${err_pos.join(', ')}]` : "No unknown errors found."
    ];
    
    return { success: true, text: parsed_text, err_pos, magnitudes, synd, bm_steps, logs };
}

// ==========================================
// 6. UI INTERACTIVE ROUTINES
// ==========================================
// Global state
let originalMsg = [];
let encodedMsg = [];
let rsEccCount = 6;

let currentQRGrid = [];
let currentQRTypes = [];
let currentQRVersion = 1;
let currentQRSize = 21;
let currentECLevel = 'M';
let currentQRMask = 0;

let simGrid = [];
let simTypes = [];
let isDrawing = false;

// Initialize Tab Buttons
document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
        document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
        document.querySelectorAll('.tab-pane').forEach(p => p.classList.remove('active'));
        btn.classList.add('active');
        const tabId = btn.getAttribute('data-tab');
        document.getElementById(tabId).classList.add('active');
    });
});

// GF(256) Multiplication
document.getElementById('btn-calc-mul').addEventListener('click', () => {
    const a = parseInt(document.getElementById('val-a').value) || 0;
    const b = parseInt(document.getElementById('val-b').value) || 0;
    
    const prod = gf_mul(a, b);
    const steps = [
        `Input A = ${a} <=> ${byte_to_poly_str(a)}`,
        `Input B = ${b} <=> ${byte_to_poly_str(b)}`,
        a === 0 || b === 0 ? "One of the inputs is 0. Product is 0." : `Log(A) = ${gf_log[a]}`,
        a === 0 || b === 0 ? "" : `Log(B) = ${gf_log[b]}`,
        a === 0 || b === 0 ? "" : `Sum of Logs = ${gf_log[a]} + ${gf_log[b]} = ${gf_log[a] + gf_log[b]}`,
        `Product = ${prod} <=> ${byte_to_poly_str(prod)}`
    ].filter(Boolean);
    
    document.getElementById('gf-mul-output').style.display = 'block';
    document.getElementById('eq-a').textContent = `A = ${a} <=> ${byte_to_poly_str(a)}`;
    document.getElementById('eq-b').textContent = `B = ${b} <=> ${byte_to_poly_str(b)}`;
    document.getElementById('eq-prod').textContent = `Result = ${prod} <=> ${byte_to_poly_str(prod)}`;
    
    const stepsList = document.getElementById('gf-mul-steps');
    stepsList.innerHTML = '';
    steps.forEach(step => {
        const li = document.createElement('li');
        li.textContent = step;
        stepsList.appendChild(li);
    });
});

// Polynomial Multiplication
document.getElementById('btn-poly-mul').addEventListener('click', () => {
    const p1 = document.getElementById('poly-a').value.trim().split(/\s+/).map(Number);
    const p2 = document.getElementById('poly-b').value.trim().split(/\s+/).map(Number);
    
    const result = poly_mul(p1, p2);
    
    document.getElementById('poly-mul-output').style.display = 'block';
    document.getElementById('poly-coeff-res').textContent = result.join(' ');
    
    const stepsList = document.getElementById('poly-mul-steps');
    stepsList.innerHTML = '';
    
    let stepLog = [
        `P1 = ` + p1.map((v, i) => `${v}*x^${p1.length - 1 - i}`).join(" + "),
        `P2 = ` + p2.map((v, i) => `${v}*x^${p2.length - 1 - i}`).join(" + "),
        "All products computed: "
    ];
    
    for (let i = 0; i < p1.length; i++) {
        for (let j = 0; j < p2.length; j++) {
            if (p1[i] !== 0 && p2[j] !== 0) {
                let termProd = gf_mul(p1[i], p2[j]);
                let pow = (p1.length - 1 - i) + (p2.length - 1 - j);
                stepLog.push(`  (${p1[i]} * ${p2[j]})x^${pow} = ${termProd}x^${pow}`);
            }
        }
    }
    stepLog.push(`Summed coefficients = [${result.join(', ')}]`);
    
    stepLog.forEach(step => {
        const li = document.createElement('li');
        li.textContent = step;
        stepsList.appendChild(li);
    });
});

// GF(256) Generator
document.getElementById('btn-construct-gf').addEventListener('click', () => {
    const prim = parseInt(document.getElementById('prim-poly-input').value) || 285;
    let exp = new Array(512).fill(0);
    let log = new Array(256).fill(0);
    let x = 1;
    let steps = [];
    
    for (let i = 0; i < 255; i++) {
        exp[i] = x;
        log[x] = i;
        steps.push(`alpha^${i} = ${x} (hex ${x.toString(16)} <=> ${byte_to_poly_str(x)})`);
        x <<= 1;
        if (x & 0x100) x ^= prim;
    }
    
    let seen = new Set();
    let duplicates = [];
    for (let val of exp.slice(0, 255)) {
        if (seen.has(val)) duplicates.push(val);
        seen.add(val);
    }
    
    const statusLbl = document.getElementById('gf-construct-status');
    document.getElementById('gf-construct-output').style.display = 'block';
    
    if (duplicates.length > 0) {
        statusLbl.textContent = `Invalid Generator. Duplicates found: ${duplicates.slice(0, 5).join(', ')}`;
        statusLbl.style.color = '#ef4444';
    } else {
        statusLbl.textContent = "Valid Galois Field GF(2^8) constructed!";
        statusLbl.style.color = '#10b981';
    }
    
    const stepsList = document.getElementById('gf-construct-steps');
    stepsList.innerHTML = '';
    steps.slice(0, 50).forEach(s => {
        const li = document.createElement('li');
        li.textContent = s;
        stepsList.appendChild(li);
    });
    if (steps.length > 50) {
        const li = document.createElement('li');
        li.textContent = `... truncated (${steps.length - 50} more steps)`;
        stepsList.appendChild(li);
    }
    
    const previewDiv = document.getElementById('gf-lut-preview');
    let previewHtml = '<strong>Log Table (gf_log):</strong><br>';
    log.forEach((val, idx) => {
        if (idx > 0 && idx < 12) previewHtml += `log[${idx}] = ${val}<br>`;
    });
    previewHtml += '...<br><br><strong>Exp Table (gf_exp):</strong><br>';
    exp.slice(0, 12).forEach((val, idx) => {
        previewHtml += `exp[${idx}] = alpha^${idx} = ${val}<br>`;
    });
    previewHtml += '...';
    previewDiv.innerHTML = previewHtml;
});

// Reed-Solomon Sandbox Encode
document.getElementById('btn-rs-encode').addEventListener('click', () => {
    const inputVal = document.getElementById('rs-input').value.trim();
    originalMsg = inputVal.split(/\s+/).map(Number);
    rsEccCount = parseInt(document.getElementById('rs-ecc-count').value) || 6;
    
    encodedMsg = rs_encode_msg(originalMsg, rsEccCount);
    
    document.getElementById('rs-encoded-section').style.display = 'block';
    renderRSBytes();
    document.getElementById('rs-decoder-output').style.display = 'none';
});

function renderRSBytes() {
    const container = document.getElementById('rs-bytes-display');
    container.innerHTML = '';
    
    encodedMsg.forEach((byte, idx) => {
        const block = document.createElement('div');
        block.className = 'byte-block';
        if (idx >= originalMsg.length) {
            block.classList.add('ecc-byte');
        }
        
        block.innerHTML = `
            <span class="val">${byte}</span>
            <span class="lbl">${idx < originalMsg.length ? 'Data' : 'ECC'}[${idx}]</span>
        `;
        
        block.addEventListener('click', () => {
            const newVal = prompt(`Enter value for byte at index ${idx} (0-255):`, byte);
            if (newVal !== null) {
                const parsed = parseInt(newVal);
                if (!isNaN(parsed) && parsed >= 0 && parsed <= 255) {
                    encodedMsg[idx] = parsed;
                    block.querySelector('.val').textContent = parsed;
                    if (parsed !== (idx < originalMsg.length ? originalMsg[idx] : dataByteAt(idx))) {
                        block.classList.add('corrupted');
                    } else {
                        block.classList.remove('corrupted');
                    }
                }
            }
        });
        container.appendChild(block);
    });
}

function dataByteAt(idx) {
    let fresh = rs_encode_msg(originalMsg, rsEccCount);
    return fresh[idx];
}

// Inject random errors
document.getElementById('btn-inject-errors').addEventListener('click', () => {
    const maxErrors = Math.floor(rsEccCount / 2);
    let count = 0;
    const blocks = document.querySelectorAll('#rs-bytes-display .byte-block');
    
    // Reset to clean
    encodedMsg = rs_encode_msg(originalMsg, rsEccCount);
    blocks.forEach((b, i) => {
        b.querySelector('.val').textContent = encodedMsg[i];
        b.classList.remove('corrupted');
    });
    
    while (count < maxErrors) {
        const idx = Math.floor(Math.random() * encodedMsg.length);
        const randomVal = Math.floor(Math.random() * 256);
        if (encodedMsg[idx] !== randomVal) {
            encodedMsg[idx] = randomVal;
            blocks[idx].querySelector('.val').textContent = randomVal;
            blocks[idx].classList.add('corrupted');
            count++;
        }
    }
});

// Reed-Solomon Decode
document.getElementById('btn-rs-decode').addEventListener('click', () => {
    const synd = rs_calc_syndromes(encodedMsg, rsEccCount);
    const bmTableContainer = document.getElementById('rs-bm-table-container');
    const bmTableBody = document.querySelector('#rs-bm-table tbody');
    bmTableBody.innerHTML = '';
    
    const outputBox = document.getElementById('rs-decoder-output');
    outputBox.style.display = 'block';
    
    const statusDiv = document.getElementById('rs-recovery-status');
    statusDiv.innerHTML = '';
    
    try {
        let [err_loc, bm_steps] = rs_find_error_locator_verbose(synd, rsEccCount);
        
        bmTableContainer.style.display = 'block';
        bm_steps.forEach(step => {
            const tr = document.createElement('tr');
            const formatPoly = (coeffs) => {
                return coeffs.map((c, i) => c === 0 ? '' : (i === 0 ? `${c}` : `${c}x^${i}`)).filter(Boolean).join(' + ') || '0';
            };
            tr.innerHTML = `
                <td>${step.r}</td>
                <td>${step.d}</td>
                <td>${step.L}</td>
                <td><code>${formatPoly(step.lam)}</code></td>
                <td><code>${formatPoly(step.B)}</code></td>
                <td style="font-family: sans-serif; font-size: 0.8rem; color: var(--text-muted);">${step.description}</td>
            `;
            bmTableBody.appendChild(tr);
        });
        
        let [corrected, err_pos, magnitudes] = rs_correct_msg(encodedMsg, rsEccCount);
        
        statusDiv.className = 'recovery-status status-success';
        statusDiv.textContent = `Success: Recovered message: [${corrected.slice(0, originalMsg.length).join(', ')}]`;
        encodedMsg = corrected;
        renderRSBytes();
        
        const errInfo = document.getElementById('rs-err-info');
        errInfo.innerHTML = '';
        if (err_pos.length > 0) {
            err_pos.forEach((pos, idx) => {
                const li = document.createElement('li');
                li.textContent = `Error position: Byte ${pos} | Magnitude: ${magnitudes[idx]}`;
                errInfo.appendChild(li);
            });
        } else {
            errInfo.innerHTML = '<li>No errors detected.</li>';
        }
    } catch (e) {
        statusDiv.className = 'recovery-status status-failed';
        statusDiv.textContent = `Failure: ${e.message}`;
    }
    
    const syndList = document.getElementById('rs-syndromes-list');
    syndList.innerHTML = '';
    synd.forEach((s, idx) => {
        if (idx === 0) return;
        const li = document.createElement('li');
        li.textContent = `S_${idx} = R(alpha^${idx}) = ${s}`;
        syndList.appendChild(li);
    });
});

// Update capacity tip
document.getElementById('rs-ecc-count').addEventListener('input', (e) => {
    const nsym = parseInt(e.target.value) || 0;
    document.getElementById('rs-correct-capacity').textContent = Math.floor(nsym / 2);
});

// QR Code Generator
document.getElementById('btn-generate-qr').addEventListener('click', () => {
    const text = document.getElementById('qr-text').value;
    const ec_level = document.getElementById('qr-ec').value;
    
    try {
        let [matrix, mask, full_bytes, mask_details] = create_qr_code_verbose(text, ec_level);
        
        currentQRGrid = matrix.grid;
        currentQRTypes = matrix.types;
        currentQRSize = matrix.size;
        currentQRVersion = matrix.version;
        currentECLevel = ec_level;
        currentQRMask = mask;
        
        renderQRGrid('qr-grid-container', currentQRGrid, currentQRTypes);
        
        document.getElementById('lbl-qr-version').textContent = `${matrix.version} (${matrix.size}x${matrix.size})`;
        document.getElementById('lbl-qr-ecc').textContent = full_bytes.length - QR_SPECS[`${matrix.version}_${ec_level}`][0];
        document.getElementById('lbl-qr-mask').textContent = mask;
        
        document.getElementById('qr-legend-box').style.display = 'block';
        document.getElementById('qr-result-box').style.display = 'flex';
        
        // Render mask comparison table
        const maskComparisonBox = document.getElementById('qr-mask-comparison-box');
        const maskTableBody = document.querySelector('#qr-mask-table tbody');
        maskTableBody.innerHTML = '';
        maskComparisonBox.style.display = 'block';
        
        const formulas = {
            0: "(row + col) % 2 == 0",
            1: "row % 2 == 0",
            2: "col % 3 == 0",
            3: "(row + col) % 3 == 0",
            4: "(floor(row/2) + floor(col/3)) % 2 == 0",
            5: "((row * col) % 2) + ((row * col) % 3) === 0",
            6: "(((row * col) % 2) + ((row * col) % 3)) % 2 == 0",
            7: "(((row + col) % 2) + ((row * col) % 3)) % 2 == 0"
        };
        
        for (let maskIdx = 0; maskIdx < 8; maskIdx++) {
            const details = mask_details[maskIdx];
            const tr = document.createElement('tr');
            if (maskIdx === mask) tr.className = 'selected-mask-row';
            
            tr.innerHTML = `
                <td>Mask ${maskIdx}</td>
                <td style="font-size: 0.8rem; font-family: monospace;">${formulas[maskIdx]}</td>
                <td>${details.n1}</td>
                <td>${details.n2}</td>
                <td>${details.n3}</td>
                <td>${details.n4}</td>
                <td style="font-weight: 700;">${details.total}</td>
                <td>${maskIdx === mask ? '✅ Selected (Lowest)' : 'Rejected'}</td>
            `;
            maskTableBody.appendChild(tr);
        }
        
        resetSimGrid();
    } catch (e) {
        alert(e.message);
    }
});

function renderQRGrid(containerId, grid, types) {
    const container = document.getElementById(containerId);
    const size = grid.length;
    container.style.gridTemplateColumns = `repeat(${size}, 1fr)`;
    container.innerHTML = '';
    
    for (let r = 0; r < size; r++) {
        for (let c = 0; c < size; c++) {
            const cell = document.createElement('div');
            cell.className = 'qr-cell';
            if (grid[r][c] === 1) cell.classList.add('dark');
            const type = types[r][c];
            if (type && type !== 'empty') cell.classList.add(`legend-${type}`);
            container.appendChild(cell);
        }
    }
}

// Simulator Reset & Scratch
function resetSimGrid() {
    simGrid = currentQRGrid.map(row => [...row]);
    simTypes = currentQRTypes.map(row => [...row]);
    renderSimGrid();
    document.getElementById('sim-status').textContent = 'Clean / Loaded';
    document.getElementById('sim-decoded-text').textContent = '-';
    document.getElementById('sim-logs-box').style.display = 'none';
    const simBmTableContainer = document.getElementById('sim-bm-table-container');
    if (simBmTableContainer) simBmTableContainer.style.display = 'none';
}

function renderSimGrid() {
    const container = document.getElementById('sim-grid-container');
    const size = simGrid.length;
    container.style.gridTemplateColumns = `repeat(${size}, 1fr)`;
    container.innerHTML = '';
    
    for (let r = 0; r < size; r++) {
        for (let c = 0; c < size; c++) {
            const cell = document.createElement('div');
            cell.className = 'qr-cell';
            if (simGrid[r][c] === 1) cell.classList.add('dark');
            else if (simGrid[r][c] === -1) cell.classList.add('scratched');
            
            cell.addEventListener('mousedown', () => { isDrawing = true; scratchCell(r, c); });
            cell.addEventListener('mouseenter', () => { if (isDrawing) scratchCell(r, c); });
            container.appendChild(cell);
        }
    }
}
window.addEventListener('mouseup', () => { isDrawing = false; });

function scratchCell(r, c) {
    const penSize = parseInt(document.getElementById('scratch-size').value) || 2;
    const size = simGrid.length;
    for (let dr = -penSize + 1; dr < penSize; dr++) {
        for (let dc = -penSize + 1; dc < penSize; dc++) {
            const nr = r + dr;
            const nc = c + dc;
            if (nr >= 0 && nr < size && nc >= 0 && nc < size) {
                simGrid[nr][nc] = -1; // scratched
            }
        }
    }
    renderSimGrid();
}

document.getElementById('btn-sim-reset').addEventListener('click', resetSimGrid);
document.getElementById('btn-sim-scratch').addEventListener('click', () => {
    const size = simGrid.length;
    const count = Math.floor(size * size * 0.20);
    for (let i = 0; i < count; i++) {
        const r = Math.floor(Math.random() * size);
        const c = Math.floor(Math.random() * size);
        simGrid[r][c] = -1;
    }
    renderSimGrid();
});

// Simulator Decode
document.getElementById('btn-sim-decode').addEventListener('click', () => {
    const statusLbl = document.getElementById('sim-status');
    const textLbl = document.getElementById('sim-decoded-text');
    const logsBox = document.getElementById('sim-logs-box');
    const logsList = document.getElementById('sim-logs-list');
    
    logsBox.style.display = 'block';
    logsList.innerHTML = '';
    
    try {
        let result = decode_qr_grid_verbose(simGrid, simTypes);
        statusLbl.textContent = 'Successfully Repaired & Decoded!';
        statusLbl.style.color = '#10b981';
        textLbl.textContent = result.text;
        
        result.logs.forEach(log => {
            const li = document.createElement('li');
            li.textContent = log;
            logsList.appendChild(li);
        });
        
        // Render BM steps
        const simBmTableContainer = document.getElementById('sim-bm-table-container');
        const simBmTableBody = document.querySelector('#sim-bm-table tbody');
        simBmTableBody.innerHTML = '';
        simBmTableContainer.style.display = 'block';
        
        result.bm_steps.forEach(step => {
            const tr = document.createElement('tr');
            const formatPoly = (coeffs) => {
                return coeffs.map((c, i) => c === 0 ? '' : (i === 0 ? `${c}` : `${c}x^${i}`)).filter(Boolean).join(' + ') || '0';
            };
            tr.innerHTML = `
                <td>${step.r}</td>
                <td>${step.d}</td>
                <td>${step.L}</td>
                <td><code>${formatPoly(step.lam)}</code></td>
                <td><code>${formatPoly(step.B)}</code></td>
                <td style="font-family: sans-serif; font-size: 0.8rem; color: var(--text-muted);">${step.description}</td>
            `;
            simBmTableBody.appendChild(tr);
        });
    } catch (e) {
        statusLbl.textContent = 'Decoding Failed';
        statusLbl.style.color = '#ef4444';
        textLbl.textContent = '-';
        
        const li = document.createElement('li');
        li.textContent = `Error: ${e.message}`;
        logsList.appendChild(li);
        
        const li2 = document.createElement('li');
        li2.textContent = "Damage exceeded Reed-Solomon correction capacity.";
        logsList.appendChild(li2);
    }
});

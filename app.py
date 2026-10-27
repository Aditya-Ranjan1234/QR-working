# app.py
from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import Optional
import os

from gf256 import gf_mul, gf_add, gf_div, gf_inv, gf_exp, gf_log
from polynomial import poly_mul, poly_add, poly_eval
from reedsolomon import rs_encode_msg, rs_calc_syndromes, rs_find_error_locator, rs_find_error_locator_verbose, rs_find_errors, rs_forney_magnitudes, rs_correct_msg
from qr_encoder import create_qr_code, create_qr_code_verbose, QR_SPECS
from qr_decoder import decode_qr_grid

app = FastAPI(title="QR & Reed-Solomon Scratch Visualizer")

# Request Models
class GF256MulReq(BaseModel):
    a: int
    b: int

class PolyMulReq(BaseModel):
    p1: list[int]
    p2: list[int]

class QREncodeReq(BaseModel):
    text: str
    ec_level: str

class QRDecodeReq(BaseModel):
    grid: list[list[int]]
    types: list[list[str]]
    erased_codewords: Optional[list[int]] = None

class RSEncodeReq(BaseModel):
    msg: list[int]
    nsym: int

class RSDecodeReq(BaseModel):
    msg: list[int]
    nsym: int
    erasures: Optional[list[int]] = None

class GF256TableReq(BaseModel):
    prim: int

# Math Helper: format byte as polynomial string
def byte_to_poly_str(val: int) -> str:
    if val == 0:
        return "0"
    terms = []
    for i in range(7, -1, -1):
        if (val >> i) & 1:
            if i == 0:
                terms.append("1")
            elif i == 1:
                terms.append("x")
            else:
                terms.append(f"x^{i}")
    return " + ".join(terms)

@app.post("/api/math/gf256-mul")
def api_gf256_mul(req: GF256MulReq):
    a = req.a & 0xFF
    b = req.b & 0xFF
    
    if a == 0 or b == 0:
        result = 0
        steps = [f"One of the inputs is 0. Product is 0."]
    else:
        log_a = gf_log[a]
        log_b = gf_log[b]
        sum_log = log_a + log_b
        result = gf_mul(a, b)
        steps = [
            f"Input A = {a} ({bin(a)} <=> {byte_to_poly_str(a)})",
            f"Input B = {b} ({bin(b)} <=> {byte_to_poly_str(b)})",
            f"Log(A) = {log_a} (alpha^{log_a} = {a})",
            f"Log(B) = {log_b} (alpha^{log_b} = {b})",
            f"Sum of Logs = {log_a} + {log_b} = {sum_log}",
            f"Result = exp({sum_log}) = exp({sum_log % 255}) = {result} ({bin(result)} <=> {byte_to_poly_str(result)})"
        ]
        
    return {
        "a": a,
        "b": b,
        "a_poly": byte_to_poly_str(a),
        "b_poly": byte_to_poly_str(b),
        "result": result,
        "result_poly": byte_to_poly_str(result),
        "steps": steps
    }

@app.post("/api/math/poly-mul")
def api_poly_mul(req: PolyMulReq):
    p1 = req.p1
    p2 = req.p2
    result = poly_mul(p1, p2)
    
    # Generate human readable steps
    steps = []
    steps.append(f"P1 = " + " | ".join(f"{x}*x^{len(p1)-1-i}" for i, x in enumerate(p1)))
    steps.append(f"P2 = " + " | ".join(f"{x}*x^{len(p2)-1-i}" for i, x in enumerate(p2)))
    
    components = []
    for i, c1 in enumerate(p1):
        for j, c2 in enumerate(p2):
            if c1 != 0 and c2 != 0:
                prod = gf_mul(c1, c2)
                pow_idx = (len(p1) - 1 - i) + (len(p2) - 1 - j)
                components.append(f"({c1} * {c2})x^{pow_idx} = {prod}x^{pow_idx}")
                
    steps.append("All non-zero term multiplications:")
    steps.extend(components)
    steps.append(f"Result polynomial coefficients (highest degree first): {result}")
    
    return {
        "p1": p1,
        "p2": p2,
        "result": result,
        "steps": steps
    }

@app.post("/api/encode")
def api_encode_qr(req: QREncodeReq):
    try:
        matrix, mask, full_bytes, mask_details = create_qr_code_verbose(req.text, req.ec_level)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
        
    version = matrix.version
    data_cap, ecc_cap, size = QR_SPECS[(version, req.ec_level)]
    
    data_bytes = full_bytes[:data_cap]
    ecc_bytes = full_bytes[data_cap:]
    
    return {
        "version": version,
        "size": size,
        "ec_level": req.ec_level,
        "mask_pattern": mask,
        "data_bytes": data_bytes,
        "ecc_bytes": ecc_bytes,
        "grid": matrix.grid,
        "types": matrix.types,
        "mask_details": mask_details,
    }

@app.post("/api/decode")
def api_decode_qr(req: QRDecodeReq):
    try:
        erased_idx = req.erased_codewords or []
        text, err_pos, magnitudes = decode_qr_grid(req.grid, req.types, erased_idx)
        
        # Verbose trace for step-by-step UI table
        import qr_decoder
        left_bits, right_bits = qr_decoder.extract_format_bits(req.grid)
        ec_level, mask_pattern = qr_decoder.decode_format_info(left_bits, right_bits)
        version = 1 if len(req.grid) == 21 else 2
        total_data_ecc_bytes = 26 if version == 1 else 44
        data_cap, ecc_bytes_count, _ = QR_SPECS[(version, ec_level)]
        
        unmasked_bits, dynamic_erasures = qr_decoder.extract_and_unmask_data_bits(req.grid, req.types, mask_pattern)
        raw_bytes = qr_decoder.bits_to_bytes(unmasked_bits)[:total_data_ecc_bytes]
        
        synd = rs_calc_syndromes(raw_bytes, ecc_bytes_count)
        
        combined_erasures = sorted(list(set(erased_idx + dynamic_erasures)))
        combined_erasures = [e for e in combined_erasures if e < total_data_ecc_bytes]
        combined_erasures_end = [len(raw_bytes) - 1 - e for e in combined_erasures]
        
        _, bm_steps = rs_find_error_locator_verbose(synd, ecc_bytes_count, combined_erasures_end)
        
        logs = []
        logs.append(f"Detected QR Version: {version} | EC Level: {ec_level} | Mask Pattern: {mask_pattern}")
        logs.append(f"Extracted {len(raw_bytes)} codeword bytes from zig-zag grid.")
        if combined_erasures:
            logs.append(f"Scratches detected at codeword bytes: {combined_erasures} (treated as erasures).")
        if err_pos:
            logs.append(f"Reed-Solomon solved and corrected {len(err_pos)} error(s) at positions: {err_pos}")
            for pos, mag in zip(err_pos, magnitudes):
                logs.append(f" - Byte {pos} corrected with magnitude {mag}")
        else:
            logs.append("No unknown byte errors found in data.")
            
        return {
            "success": True,
            "text": text,
            "error_positions": err_pos,
            "magnitudes": magnitudes,
            "syndromes": synd,
            "bm_steps": bm_steps,
            "logs": logs
        }
    except Exception as e:
        # Fallback trace in case of complete decode failure
        try:
            import qr_decoder
            left_bits, right_bits = qr_decoder.extract_format_bits(req.grid)
            ec_level, mask_pattern = qr_decoder.decode_format_info(left_bits, right_bits)
            version = 1 if len(req.grid) == 21 else 2
            total_data_ecc_bytes = 26 if version == 1 else 44
            data_cap, ecc_bytes_count, _ = QR_SPECS[(version, ec_level)]
            unmasked_bits, dynamic_erasures = qr_decoder.extract_and_unmask_data_bits(req.grid, req.types, mask_pattern)
            raw_bytes = qr_decoder.bits_to_bytes(unmasked_bits)[:total_data_ecc_bytes]
            synd = rs_calc_syndromes(raw_bytes, ecc_bytes_count)
            dynamic_erasures_end = [len(raw_bytes) - 1 - e for e in dynamic_erasures]
            _, bm_steps = rs_find_error_locator_verbose(synd, ecc_bytes_count, dynamic_erasures_end)
        except Exception:
            synd = []
            bm_steps = []
        return {
            "success": False,
            "error": str(e),
            "syndromes": synd,
            "bm_steps": bm_steps,
            "logs": [f"Decoding failed: {str(e)}", "Damage exceeded Reed-Solomon correction capacity."]
        }

@app.post("/api/rs/encode")
def api_rs_encode(req: RSEncodeReq):
    try:
        encoded = rs_encode_msg(req.msg, req.nsym)
        return {"encoded": encoded}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/api/rs/decode")
def api_rs_decode(req: RSDecodeReq):
    synd = rs_calc_syndromes(req.msg, req.nsym)
    erase_loc = req.erasures or []
    erase_loc_end = [len(req.msg) - 1 - e for e in erase_loc]
    try:
        # Trace BM
        err_loc, bm_steps = rs_find_error_locator_verbose(synd, req.nsym, erase_loc_end)
        corrected, err_pos, magnitudes = rs_correct_msg(req.msg, req.nsym, erase_loc)
        return {
            "success": True,
            "corrected": corrected,
            "error_positions": err_pos,
            "magnitudes": magnitudes,
            "syndromes": synd,
            "bm_steps": bm_steps
        }
    except Exception as e:
        # Trace BM even if decoding fails downstream
        try:
            _, bm_steps = rs_find_error_locator_verbose(synd, req.nsym, erase_loc_end)
        except Exception:
            bm_steps = []
        return {
            "success": False,
            "error": str(e),
            "syndromes": synd,
            "error_positions": [],
            "magnitudes": [],
            "bm_steps": bm_steps
        }

@app.post("/api/math/gf256-table")
def api_gf256_table(req: GF256TableReq):
    prim = req.prim
    exp = [0] * 512
    log = [0] * 256
    x = 1
    conflict = False
    history = []
    
    for i in range(255):
        exp[i] = x
        log[x] = i
        history.append(f"alpha^{i} = {x} (hex {hex(x)} <=> {byte_to_poly_str(x)})")
        x <<= 1
        if x & 0x100:
            x ^= prim
            
    seen = set()
    duplicates = []
    for val in exp[:255]:
        if val in seen:
            duplicates.append(val)
        seen.add(val)
        
    if duplicates:
        conflict = True
        status = f"Invalid Primitive Polynomial. Duplicates: {duplicates[:5]}... Generated field contains {len(seen)} unique elements."
    else:
        status = "Valid Galois Field GF(2^8) constructed!"
        
    return {
        "status": status,
        "is_valid": not conflict,
        "history": history[:50] + ([f"... truncated ({len(history) - 50} more steps)"] if len(history) > 50 else []),
        "exp": exp[:255],
        "log": log
    }

# Serve Frontend static files
os.makedirs("static", exist_ok=True)
app.mount("/", StaticFiles(directory="static", html=True), name="static")

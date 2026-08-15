# ============================================================================
# Module 08 - Drug target: mini-tabix utility
# Fetch a genomic region from a remote bgzf+tabix file via HTTP range
# requests (used to pull eQTL Catalogue GTEx v8 per-tissue cis summary stats
# for the SERPINE1 locus).
# Usage: py -3 t1_minitabix.py <base_url_without_extension> <chrom> <beg> <end> <out.tsv>
# ============================================================================

import os as _os
# --- repository paths -----------------------------------------------------------
# REPO resolves to the repository root (code/<module>/this_file.py -> ../..).
# Override with the REPRO_ROOT environment variable if you run from elsewhere.
REPO = _os.environ.get("REPRO_ROOT", _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..")))
# GBD_RAW_DIR: folder with the raw IHME GBD 2021 csv downloads (see data/README.md)
GBD_RAW_DIR = _os.environ.get("GBD_RAW_DIR", _os.path.join(REPO, "data", "raw", "gbd_2021"))
#!/usr/bin/env python3
"""Mini-tabix: fetch a genomic region from a remote bgzf+tabix file via HTTP range requests.
Usage: py -3 t1_minitabix.py <base_url_without_extension> <chrom> <beg> <end> <out.tsv>
e.g.:  py -3 t1_minitabix.py https://.../QTD000271/QTD000271.all.tsv.gz 7 100126167 101138431 out.tsv
"""
import sys, gzip, struct, zlib, urllib.request, subprocess, tempfile, os

def fetch_range(url, start, end):
    for attempt in range(3):
        try:
            with tempfile.NamedTemporaryFile(delete=False) as tf:
                tmp = tf.name
            subprocess.run(["curl", "-s", "-m", "180", "-o", tmp, "-r", f"{start}-{end}", url],
                           check=True, timeout=200)
            with open(tmp, "rb") as fh:
                data = fh.read()
            os.unlink(tmp)
            return data
        except Exception as e:
            print(f"range fetch attempt {attempt+1} failed: {e}", flush=True)
    raise RuntimeError("range fetch failed")

def fetch_full(url):
    with tempfile.NamedTemporaryFile(delete=False) as tf:
        tmp = tf.name
    subprocess.run(["curl", "-s", "-m", "180", "-o", tmp, url], check=True, timeout=200)
    with open(tmp, "rb") as fh:
        data = fh.read()
    os.unlink(tmp)
    return data

def parse_tbi(buf):
    assert buf[:4] == b"TBI\x01", "not a tbi"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", buf[4:36])
    names = buf[36:36 + l_nm].split(b"\x00")[:-1]
    off = 36 + l_nm
    refs = {}
    for r in range(n_ref):
        n_bin = struct.unpack("<i", buf[off:off+4])[0]; off += 4
        bins = {}
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack("<Ii", buf[off:off+8]); off += 8
            chunks = []
            for _ in range(n_chunk):
                beg, end = struct.unpack("<QQ", buf[off:off+16]); off += 16
                chunks.append((beg, end))
            bins[bin_id] = chunks
        n_intv = struct.unpack("<i", buf[off:off+4])[0]; off += 4 + 8 * n_intv
        refs[names[r].decode()] = bins
    return refs

def reg2bins(beg, end):
    end -= 1
    bins = [0]
    for k in range(1 + (beg >> 26), 2 + (end >> 26)): bins.append(k)
    for k in range(9 + (beg >> 23), 10 + (end >> 23)): bins.append(k)
    for k in range(73 + (beg >> 20), 74 + (end >> 20)): bins.append(k)
    for k in range(585 + (beg >> 17), 586 + (end >> 17)): bins.append(k)
    for k in range(4681 + (beg >> 14), 4682 + (end >> 14)): bins.append(k)
    return bins

def main():
    url, chrom, beg, end, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    tbi_raw = fetch_full(url + ".tbi")
    refs = parse_tbi(gzip.decompress(tbi_raw))
    key = chrom if chrom in refs else ("chr" + chrom if "chr" + chrom in refs else None)
    if key is None:
        print("chrom not found; available:", list(refs)[:5]); sys.exit(1)
    chunks = []
    for b in reg2bins(beg, end):
        chunks += refs[key].get(b, [])
    if not chunks:
        print("no chunks for region"); sys.exit(1)
    # byte ranges from virtual offsets; merge
    ranges = sorted((c[0] >> 16, (c[1] >> 16) + 65536) for c in chunks)
    merged = [list(ranges[0])]
    for s, e in ranges[1:]:
        if s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
        else: merged.append([s, e])
    blob = b""
    for s, e in merged:
        blob += fetch_range(url, s, e - 1)
    # decompress BGZF members individually; tolerate trailing partial block
    parts = []
    while blob:
        d = zlib.decompressobj(31)
        try:
            out_b = d.decompress(blob)
        except Exception:
            break
        parts.append(out_b)
        unused = d.unused_data
        if not unused:
            break
        blob = unused
    data = b"".join(parts)
    # header is at the start of the FILE (not in our range): fetch separately
    head_blob = fetch_range(url, 0, 65535)
    htxt = zlib.decompressobj(31).decompress(head_blob).decode("utf-8", "replace")
    header = [l for l in htxt.splitlines() if l.startswith("#") or "molecular_trait_id" in l or "chromosome" in l][-1]
    lines = data.decode("utf-8", "replace").splitlines()
    cols = header.lstrip("#").split("\t")
    try:
        ci = cols.index("chromosome")
    except ValueError:
        ci = 0
    pi = cols.index("position") if "position" in cols else (ci + 1)
    out_lines = [header]
    for l in lines[1:]:   # first line of blob may be partial
        if not l.strip(): continue
        f = l.split("\t")
        if len(f) <= max(ci, pi): continue
        try:
            pos = int(f[pi])
        except ValueError:
            continue
        if f[ci] in (chrom, "chr" + chrom) and beg <= pos <= end:
            out_lines.append(l)
    with open(out, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    print(f"{out}: {len(out_lines)-1} data rows (from {len(merged)} range(s), {len(blob)/1e6:.1f} MB fetched)")

main()

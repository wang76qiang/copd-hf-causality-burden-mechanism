#!/usr/bin/env python3
"""Extract the 10 COPD instrument SNPs from GWAS Catalog harmonised BBJ 2025 HF GWAS
(GCST90668009 all-cause HF, GCST90668010 HFrEF, GCST90668011 HFpEF) via HTTP range tabix.
"""
import json, gzip, struct, zlib, urllib.request, subprocess, tempfile, os, sys, csv

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90668001-GCST90669000"
ACCS = {
    "GCST90668009": "BBJ_allHF",
    "GCST90668010": "BBJ_HFrEF",
    "GCST90668011": "BBJ_HFpEF",
}

def fetch_range(url, start, end):
    for attempt in range(3):
        try:
            with tempfile.NamedTemporaryFile(delete=False) as tf:
                tmp = tf.name
            subprocess.run(["curl", "-s", "-m", "300", "-o", tmp, "-r", f"{start}-{end}", url],
                           check=True, timeout=320)
            with open(tmp, "rb") as fh:
                data = fh.read()
            os.unlink(tmp)
            return data
        except Exception as e:
            print(f"  range fetch attempt {attempt+1} failed: {e}", flush=True)
    raise RuntimeError("range fetch failed")

def parse_tbi(buf):
    assert buf[:4] == b"TBI\x01"
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

def bgzf_decompress(blob):
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
    return b"".join(parts)

def get_header(url):
    blob = fetch_range(url, 0, 131071)
    txt = bgzf_decompress(blob).decode("utf-8", "replace")
    for l in txt.splitlines():
        if "hm_rsid" in l or "base_pair_location" in l:
            return l.split("\t")
    raise RuntimeError("header not found")

def main():
    snps = json.load(open("snps10_hg38.json"))
    wanted = set(snps.keys())
    os.makedirs("bbj2025_extracts", exist_ok=True)
    for acc, tag in ACCS.items():
        url = f"{BASE}/{acc}/harmonised/{acc}.h.tsv.gz"
        tbi = parse_tbi(gzip.decompress(fetch_range(url + ".tbi", 0, 10**9)))
        cols = get_header(url)
        print(acc, "cols:", cols[:8], "...", flush=True)
        # locate columns
        def col(name, *alts):
            for n in (name,) + alts:
                if n in cols: return cols.index(n)
            return None
        i_rsid = col("hm_rsid", "rsid", "variant_id")
        i_chrom = col("hm_chrom", "chromosome")
        i_pos = col("hm_pos", "base_pair_location")
        rows_out = []
        for rs, (chrom, pos, _) in snps.items():
            chrom, pos = str(chrom), int(pos)
            key = chrom if chrom in tbi else ("chr" + chrom if "chr" + chrom in tbi else None)
            if key is None:
                print(f"  {rs}: chrom {chrom} not in tbi"); continue
            chunks = []
            for b in reg2bins(pos - 2000, pos + 2000):
                chunks += tbi[key].get(b, [])
            if not chunks:
                print(f"  {rs}: no chunks"); continue
            ranges = sorted((c[0] >> 16, (c[1] >> 16) + 65536) for c in chunks)
            merged = [list(ranges[0])]
            for s, e in ranges[1:]:
                if s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
                else: merged.append([s, e])
            blob = b""
            for s, e in merged:
                blob += fetch_range(url, s, e - 1)
            data = bgzf_decompress(blob).decode("utf-8", "replace")
            found = 0
            for l in data.splitlines():
                f = l.split("\t")
                if len(f) <= max(i_rsid or 0, i_chrom, i_pos): continue
                if f[i_rsid] == rs:
                    rows_out.append(f); found += 1
            if not found:
                print(f"  {rs}: NOT FOUND in {acc}")
        with open(f"bbj2025_extracts/{tag}.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(cols)
            for f in rows_out: w.writerow(f)
        print(f"{acc} ({tag}): {len(rows_out)} SNPs extracted", flush=True)

if __name__ == "__main__":
    main()

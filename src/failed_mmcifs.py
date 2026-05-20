#!/usr/bin/env python3
import argparse, os, time, requests

CIF_URL_TMPL = "https://files.rcsb.org/download/{pdb_id}.cif"

def load_ids(path):
    ids = []
    with open(path, "r") as f:
        for line in f:
            s = line.strip().upper()
            if s and s[0].isalnum():
                ids.append(s)
    return ids

def fetch_mmcif_bytes(pdb_id, session, timeout=60):
    url = CIF_URL_TMPL.format(pdb_id=pdb_id.lower())
    try:
        r = session.get(url, timeout=timeout)
        if r.status_code == 200 and r.content and len(r.content) > 1024:
            return r.content
    except requests.RequestException:
        pass
    return None

def main():
    ap = argparse.ArgumentParser(description="Retry downloading only failed PDB mmCIF IDs.")
    ap.add_argument("--in", dest="infile", required=True, help="path to text file with failed IDs (one per line)")
    ap.add_argument("--out", default="pdb_subset", help="output folder (same as your main run)")
    ap.add_argument("--retries", type=int, default=3, help="retries per ID (default 3)")
    ap.add_argument("--timeout", type=int, default=60, help="HTTP timeout (s)")
    ap.add_argument("--delay", type=float, default=0.1, help="delay between attempts (s)")
    ap.add_argument("--no-skip", action="store_true", help="force re-download even if file exists")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    ids = load_ids(args.infile)
    if not ids:
        print("No IDs found in input file.")
        return

    print(f"Retrying {len(ids)} failed ID(s) → {os.path.abspath(args.out)}")
    ok_cnt = 0
    skip_cnt = 0
    fail_cnt = 0

    with requests.Session() as s:
        s.headers.update({"User-Agent": "mmCIF-retry/1.0 (contact: your_email@example.com)"})
        for i, pdb_id in enumerate(ids, 1):
            out_path = os.path.join(args.out, f"{pdb_id.lower()}.mmcif")
            if not args.no_skip and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
                print(f"[{i}/{len(ids)}] Skip (exists): {pdb_id.lower()}.mmcif")
                skip_cnt += 1
                time.sleep(max(0.0, args.delay))
                continue

            data = None
            for attempt in range(args.retries + 1):
                data = fetch_mmcif_bytes(pdb_id, s, timeout=args.timeout)
                if data:
                    break
                if attempt < args.retries:
                    time.sleep(1.0 + 0.5 * attempt)  # backoff

            if data:
                try:
                    with open(out_path, "wb") as fh:
                        fh.write(data)
                    print(f"[{i}/{len(ids)}] Downloaded {pdb_id.lower()}.mmcif")
                    ok_cnt += 1
                except Exception as e:
                    print(f"[{i}/{len(ids)}] Write failed {pdb_id}: {e}")
                    fail_cnt += 1
            else:
                print(f"[{i}/{len(ids)}] Failed again: {pdb_id}")
                fail_cnt += 1

            time.sleep(max(0.0, args.delay))

    print("\nRetry complete.")
    print(f"  Downloaded: {ok_cnt}")
    print(f"  Skipped:    {skip_cnt}")
    print(f"  Failed:     {fail_cnt}")

if __name__ == "__main__":
    main()

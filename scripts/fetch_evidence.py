#!/usr/bin/env python3
"""Fetch public evidence for one ClinVar variant: ClinVar, gnomAD, Ensembl VEP, PubMed.

Stdlib only. Every source degrades to None (never 0) when it has no data,
because 0 is a real score and NA is not evidence either way.

Usage (all identifiers point at the SAME variant):
  python3 fetch_evidence.py \
      --clinvar 3634421 \
      --gnomad 17-7669613-T-A \
      --hgvs "NM_000546.6:c.1178A>T" --transcript ENST00000269305 \
      --pubmed-terms "p.Asp393Val,D393V" \
      --out D393V.json

Run `--selftest` to check all four endpoints answer before trusting a batch.
"""
import argparse, json, random, re, sys, time, urllib.error, urllib.parse, urllib.request
from datetime import datetime, timezone

UA = {"User-Agent": "hermes-athena-vus/1.0 (variant curation research)"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
VEP_REST = "https://rest.ensembl.org"
VEP_TRIES = 6   # Ensembl VEP returns 500/503 on ~half of requests; retry hard
TIMEOUT = 45


def _get(url, accept=None, tries=3):
    hdrs = dict(UA)
    if accept:
        hdrs["Accept"] = accept
    err = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=hdrs)
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            err = e
            # 4xx other than 429 is a real answer (bad id, bad param) - do not burn retries
            if 400 <= e.code < 500 and e.code != 429:
                raise
        except Exception as e:
            err = e
        # Ensembl VEP flaps: identical URLs return 200/500/503 within seconds.
        backoff = min(2.0 * (i + 1), 8.0) + random.uniform(0, 1.0)
        time.sleep(backoff)
    raise RuntimeError(f"GET {url} failed after {tries} tries: {err!r}")


def _post_json(url, payload, tries=3):
    err = None
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, data=json.dumps(payload).encode(),
                headers={**UA, "Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            err = e
            body = ""
            try:
                body = e.read().decode()[:300]  # GraphQL errors explain schema drift
            except Exception:
                pass
            if body:
                err = RuntimeError(f"{e!r} {body}")
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"POST {url} failed: {err!r}")


# ---------------------------------------------------------------- ClinVar
def fetch_clinvar(variation_id):
    url = f"{EUTILS}/esummary.fcgi?db=clinvar&retmode=json&id={urllib.parse.quote(str(variation_id))}"
    rec = json.loads(_get(url))["result"].get(str(variation_id), {})
    if not rec:
        return {"variation_id": variation_id, "error": "not found", "classification": None}
    gc = rec.get("germline_classification", {}) or {}
    return {
        "variation_id": str(variation_id),
        "accession": rec.get("accession"),
        "title": rec.get("title"),
        "gene_sort": rec.get("gene_sort"),
        "classification": gc.get("description"),
        "review_status": gc.get("review_status"),
        "last_evaluated": gc.get("last_evaluated"),
        "source": url,
    }


# ----------------------------------------------------------------- gnomAD
GNOMAD_Q = '''query($id:String!,$ds:DatasetId!){
  variant(variantId:$id, dataset:$ds){ variant_id rsids
    exome{ac an homozygote_count faf95{popmax popmax_population}}
    genome{ac an homozygote_count faf95{popmax popmax_population}}
    joint{ac an homozygote_count faf95{popmax popmax_population}}
  }}'''


def fetch_gnomad(variant_id, dataset="gnomad_r4"):
    if not variant_id:
        return {"variant_id": None, "note": "no gnomAD id supplied"}
    d = _post_json(GNOMAD_API, {"query": GNOMAD_Q, "variables": {"id": variant_id, "ds": dataset}})
    v = (d.get("data") or {}).get("variant")
    if not v:
        return {"variant_id": variant_id, "dataset": dataset,
                "absent": True, "af": None, "source": GNOMAD_API,
                "note": "variant not found in dataset -> candidate for PM2, but confirm the id is right"}
    out = {"variant_id": v["variant_id"], "rsids": v.get("rsids") or [],
           "dataset": dataset, "source": GNOMAD_API, "absent": False}
    for pop in ("joint", "exome", "genome"):
        blk = v.get(pop)
        if not blk:
            continue
        faf = blk.get("faf95") or {}
        ac, an = blk.get("ac"), blk.get("an")
        out[pop] = {"ac": ac, "an": an,
                    "af": (ac / an) if (ac is not None and an) else None,
                    "homozygotes": blk.get("homozygote_count"),
                    "faf95_popmax": faf.get("popmax"),
                    "popmax_population": faf.get("popmax_population")}
    primary = out.get("joint") or out.get("exome") or out.get("genome") or {}
    out["af"] = primary.get("af")
    out["homozygotes"] = primary.get("homozygotes")
    out["faf95_popmax"] = primary.get("faf95_popmax")
    out["popmax_population"] = primary.get("popmax_population")
    return out


# -------------------------------------------------------------------- VEP
# Flags are opt-in and camelCase. Lowercase flags are silently ignored by Ensembl.
VEP_SCORE_FLAGS = "CADD=1&REVEL=1&PolyPhen=1&SIFT=1&AlphaMissense=1"
VEP_DBNSFP = "dbNSFP=BayesDel_addAF_score"   # dbNSFP and AlphaMissense do not return together


def fetch_vep(hgvs, transcript=None):
    if not hgvs:
        return {"hgvs": None, "note": "no HGVS supplied"}
    # NOTE: dbNSFP is deliberately NOT on this URL. It suppresses AlphaMissense,
    # so it gets its own call below and the two results are merged.
    if transcript:
        url = f"{VEP_REST}/vep/human/hgvs/{urllib.parse.quote(hgvs, safe='')}" \
              f"?transcript_id={transcript}&{VEP_SCORE_FLAGS}&content-type=application/json"
    else:
        url = f"{VEP_REST}/vep/human/hgvs/{urllib.parse.quote(hgvs, safe='')}" \
              f"?canonical=1&{VEP_SCORE_FLAGS}&content-type=application/json"
    def _one(u):
        data = json.loads(_get(u, accept="application/json", tries=VEP_TRIES))
        return data[0] if data else {}

    d = _one(url)
    if not d:
        return {"hgvs": hgvs, "error": "VEP returned nothing", "scores": {}}
    tcs = d.get("transcript_consequences") or []
    tc = dict(tcs[0]) if tcs else {}
    # BayesDel is not a VEP plugin: it needs the dbNSFP param, which suppresses
    # AlphaMissense on the same request. Second call, merged in.
    # A failed merge must never look like "this position has no dbNSFP data": count3
    # counts only tools that returned a value, so a silently missing BayesDel shrinks
    # the denominator and can flip PP3/BP4. Record the outcome instead of guessing.
    bayesdel_status = "not requested (no transcript)" if not transcript else "not returned"
    if transcript:
        u2 = (f"{VEP_REST}/vep/human/hgvs/{urllib.parse.quote(hgvs, safe='')}"
              f"?transcript_id={transcript}&{VEP_DBNSFP}&content-type=application/json")
        try:
            d2 = _one(u2)
            t2 = (d2.get("transcript_consequences") or [{}])[0]
            if t2.get("bayesdel_addaf_score") is not None:
                tc["bayesdel_addaf_score"] = t2["bayesdel_addaf_score"]
                bayesdel_status = "ok"
        except Exception as e:                       # endpoint flapped past the retries
            bayesdel_status = f"merge failed: {str(e)[:120]}"
    scores = {
        "revel": tc.get("revel"),
        "cadd_phred": tc.get("cadd_phred"),
        "bayesdel_addaf": tc.get("bayesdel_addaf_score"),
        "polyphen": tc.get("polyphen_prediction"),   # absent == no coverage at this position
        "sift": tc.get("sift_prediction"),
        "alphamissense": (tc.get("alphamissense") or {}).get("am_class")
        if isinstance(tc.get("alphamissense"), dict) else tc.get("alphamissense"),
    }
    return {"hgvs": hgvs, "transcript_used": tc.get("transcript_id"),
            "most_severe_consequence": d.get("most_severe_consequence"),
            # carried so check_identity() can prove gnomAD and VEP are the same locus
            "seq_region_name": d.get("seq_region_name"), "start": d.get("start"),
            "allele_string": d.get("allele_string"), "strand": d.get("strand"),
            "hgvsp": d.get("hgvsp"), "hgvsc": d.get("hgvsc"),
            "scores": scores, "bayesdel_merge": bayesdel_status,
            "tools_with_data": sorted(k for k, v in scores.items() if v is not None),
            "tools_missing": sorted(k for k, v in scores.items() if v is None),
            "source": url}


# ----------------------------------------------------------------- PubMed
def fetch_pubmed(terms, retmax=20):
    """Tight phrase search first. Tight query returning 0 is a RESULT, not a failure:
    it means no publication uses that protein change, so PS3/BS3 cannot be assigned."""
    if not terms:
        return {"terms": [], "note": "no PubMed terms supplied"}
    termlist = [t.strip() for t in terms.split(",") if t.strip()] or [terms]
    q = " OR ".join(f'"{t}"[tiab]' for t in termlist)
    url = f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax={retmax}&term={urllib.parse.quote(q)}"
    r = json.loads(_get(url))["esearchresult"]
    out = {"terms": termlist, "query": q, "count": int(r.get("count", 0)),
           "pmids": r.get("idlist") or [], "source": url}
    if out["pmids"]:
        ids = ",".join(out["pmids"][:retmax])
        s = json.loads(_get(f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id={ids}"))["result"]
        out["hits"] = [{"pmid": i, "date": s.get(i, {}).get("pubdate"),
                        "title": s.get(i, {}).get("title")} for i in out["pmids"]]
    out["action_required"] = (
        "Read abstracts by hand. PS3/BS3 may ONLY be assigned from a functional assay "
        "reported in a paper you actually read. A title is not evidence."
        if out["count"] else
        "Zero hits for the exact protein change: no PS3/BS3. Try a looser gene+domain "
        "query and read by eye before concluding nothing exists.")
    return out


# ------------------------------------------------- identifier consistency guard
_C_HGVS = re.compile(r"(c\.[0-9_]+[ACGTacgt]*>[ACGTacgt]+|c\.[0-9_]+del[ACGT]*|c\.[0-9_]+ins[ACGT]+)")
_P_HGVS = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)")


def check_identity(clinvar, gnomad, vep, hgvs=None):
    """The three identifiers must describe ONE variant.

    Nothing downstream can detect a mismatch: ClinVar will happily report a
    pathogenic classification for variant A while gnomAD and VEP answer for
    variant B, and every code still fires. Refuse to blend them.
    """
    problems, seen = [], {}
    if hgvs:
        m = _C_HGVS.search(hgvs)
        seen["requested_hgvs"] = m.group(1) if m else hgvs

    if isinstance(clinvar, dict) and clinvar.get("title"):
        t = clinvar["title"]
        mc = _C_HGVS.search(t)
        seen["clinvar_title"] = t
        seen["clinvar_hgvs"] = mc.group(1) if mc else None
        if hgvs and mc and mc.group(1).lower() != seen["requested_hgvs"].lower():
            problems.append(f"ClinVar record is {mc.group(1)}, you asked for "
                            f"{seen['requested_hgvs']} -> WRONG VARIATION ID")

    if isinstance(gnomad, dict):
        gid = gnomad.get("variant_id") or ""
        seen["gnomad_id"] = gid
        gm = re.match(r"^(\w+)-(\d+)-([ACGT]+)-([ACGT]+)$", gid or "")
        if gm and isinstance(vep, dict) and vep.get("start"):
            gchr, gpos, gref, galt = gm.groups()
            vchr = str(vep.get("seq_region_name") or "")
            vpos, valleles = str(vep.get("start")), (vep.get("allele_string") or "")
            seen["vep_locus"] = f"{vchr}:{vpos} {valleles}"
            if gchr != vchr or gpos != vpos:
                problems.append(f"gnomAD id is at {gchr}:{gpos} but VEP resolved "
                                f"{vchr}:{vpos} -> DIFFERENT POSITION")
            elif "/" in valleles:
                vref, _, valt = valleles.partition("/")
                # VEP returns allele_string in transcript orientation. On a minus-strand
                # gene (TP53) that is the reverse complement of a gnomAD id, so compare
                # in genomic orientation or every minus-strand variant looks mismatched.
                if vep.get("strand") == -1:
                    _comp = str.maketrans("ACGT", "TGCA")
                    vref, valt = vref.translate(_comp), valt.translate(_comp)
                if (gref, galt) != (vref, valt):
                    problems.append(f"gnomAD alleles {gref}>{galt} vs VEP "
                                    f"{vref}>{valt} at the same position -> DIFFERENT ALT")
    return {"ok": not problems, "problems": problems, "identifiers": seen}


def selftest():
    print("ClinVar: ", bool(fetch_clinvar(3634421).get("classification")))
    print("gnomAD:  ", fetch_gnomad("17-7669613-T-A").get("af"))
    print("VEP:     ", fetch_vep("NM_000546.6:c.1178A>T", "ENST00000269305")["scores"]["revel"])
    print("PubMed:  ", fetch_pubmed("p.Asp393Val,D393V")["count"], "hits")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinvar", help="ClinVar variation id (e.g. 3634421)")
    ap.add_argument("--gnomad", help="gnomAD variant id, GRCh38, e.g. 17-7669613-T-A")
    ap.add_argument("--hgvs", help="HGVS c. notation with accession, e.g. NM_000546.6:c.1178A>T")
    ap.add_argument("--transcript", help="MANE Select ENST id (filters VEP output)")
    ap.add_argument("--pubmed-terms", help="comma-separated protein-change spellings")
    ap.add_argument("--out", help="JSON output path (default stdout)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()

    result = {"fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    for name, fn in (("clinvar", lambda: fetch_clinvar(a.clinvar)),
                     ("gnomad", lambda: fetch_gnomad(a.gnomad)),
                     ("vep", lambda: fetch_vep(a.hgvs, a.transcript)),
                     ("pubmed", lambda: fetch_pubmed(a.pubmed_terms))):
        if (name == "clinvar" and not a.clinvar) or (name == "gnomad" and not a.gnomad) \
           or (name == "vep" and not a.hgvs) or (name == "pubmed" and not a.pubmed_terms):
            result[name] = {"skipped": "identifier not supplied"}
            continue
        try:
            result[name] = fn()
        except Exception as e:                      # a dead endpoint is NA, never a zero
            result[name] = {"error": str(e)[:400]}

    result["consistency"] = check_identity(
        result.get("clinvar"), result.get("gnomad"), result.get("vep"), a.hgvs)
    if not result["consistency"]["ok"]:
        print("!! IDENTIFIER MISMATCH - DO NOT CURATE THIS PAYLOAD", file=sys.stderr)
        for pr in result["consistency"]["problems"]:
            print("!!   " + pr, file=sys.stderr)

    txt = json.dumps(result, indent=2, ensure_ascii=False)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            f.write(txt + "\n")
        print(f"wrote {a.out}")
    else:
        print(txt)


if __name__ == "__main__":
    main()

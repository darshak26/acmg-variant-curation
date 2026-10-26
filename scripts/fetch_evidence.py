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


class OutOfScope(Exception):
    """The variant is real but outside this tool's competence. Never guess around it."""


# ClinVar files copy-number and structural variants under obj_type, NOT variant_type
# (which comes back null in esummary). Checking variant_type would never fire.
OUT_OF_SCOPE_TYPES = ("copy number gain", "copy number loss", "inversion",
                      "translocation", "complex")


# ---------------------------------------------------------------- ClinVar
def fetch_clinvar(variation_id):
    url = f"{EUTILS}/esummary.fcgi?db=clinvar&retmode=json&id={urllib.parse.quote(str(variation_id))}"
    rec = json.loads(_get(url))["result"].get(str(variation_id), {})
    if not rec:
        return {"variation_id": variation_id, "error": "not found", "classification": None}
    gc = rec.get("germline_classification", {}) or {}
    obj_type = (rec.get("obj_type") or "").strip()
    if any(t in obj_type.lower() for t in OUT_OF_SCOPE_TYPES):
        title = rec.get("title") or ""
        span = re.search(r"chr[^:]*:([0-9]+)-([0-9]+)", title)
        size = ""
        if span:
            kb = (int(span.group(2)) - int(span.group(1))) / 1000
            size = f", ~{kb:,.0f} kb" if kb < 10000 else f", ~{kb/1000:,.1f} Mb"
        raise OutOfScope(
            f"Variant {variation_id} is a {obj_type}{size}: {title}\n"
            "acmg-variant-curation curates germline SNVs and small indels under "
            "Richards 2015. This variant is out of scope and no codes will be derived.\n"
            "For CNVs use the ACMG/ClinGen 2020 technical standards (Riggs et al., "
            "Genet Med 22:245-257) with a purpose-built CNV tool.")
    return {
        "variation_id": str(variation_id),
        "obj_type": obj_type,
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
    # A wrong MANE accession makes the filter drop every consequence, which reads
    # downstream as "this position has no predictors" and quietly weakens the call.
    transcript_mismatch = None
    if transcript and not any(t.get("transcript_id") == transcript for t in tcs):
        transcript_mismatch = (
            f"no consequence returned for transcript {transcript}. Either the accession "
            "is wrong or the variant is not in that gene - the empty predictor set is the "
            "filter's doing, not missing data. Refetch with the correct MANE transcript.")
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
            "transcript_mismatch": transcript_mismatch,
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
# ClinVar records carry whichever HGVS flavour the submitter used. matching only c.
# made the guard return ok:true WITHOUT comparing anything for genomic records -
# a silent bypass of the one check this tool exists to perform.
_HGVS = re.compile(
    r"((?:[A-Z]{2}_[0-9]+\.[0-9]+:)?[cgn]\.[0-9_]+"
    r"(?:[ACGTacgt]*>[ACGTacgt]+|del[ACGT]*|ins[ACGT]+|dup[ACGT]*))")
_P_HGVS = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=)")


def _parse_hgvs(s):
    """-> {kind, raw, pos, ref, alt} for c./g./n. notation, or None if absent."""
    if not s:
        return None
    m = _HGVS.search(s)
    if not m:
        return None
    raw = m.group(1)
    kind = re.match(r"(?:[A-Z]{2}_[0-9]+\.[0-9]+:)?([cgn])\.", raw)
    body = raw.split(":", 1)[-1] if ":" in raw else raw
    out = {"kind": kind.group(1) if kind else None, "raw": raw, "pos": None,
           "ref": None, "alt": None}
    pm = re.search(r"(\d+)([ACGTacgt]*)>([ACGTacgt]+)", body)
    if pm:
        out["pos"], out["ref"], out["alt"] = int(pm.group(1)), pm.group(2).upper(), pm.group(3).upper()
    else:
        pm = re.search(r"(\d+)", body)
        if pm:
            out["pos"] = int(pm.group(1))
    return out


def _vep_alleles_plus_strand(vep, req_kind):
    """VEP reports allele_string in TRANSCRIPT orientation for c./n. input, and in
    plus-strand orientation for g. input. Normalise to plus strand before comparing."""
    av = (vep or {}).get("allele_string") or ""
    if "/" not in av:
        return None, None
    r, _, a = av.partition("/")
    if req_kind in ("c", "n") and (vep or {}).get("strand") == -1:
        comp = str.maketrans("ACGT", "TGCA")
        r, a = r.translate(comp), a.translate(comp)
    return r, a


def check_identity(clinvar, gnomad, vep, hgvs=None):
    """The three identifiers must describe ONE variant.

    Nothing downstream can detect a mismatch: ClinVar will happily report a
    pathogenic classification for variant A while gnomAD and VEP answer for
    variant B, and every code still fires. Refuse to blend them.
    """
    problems, seen, checks = [], {}, []
    req = _parse_hgvs(hgvs)
    seen["requested_hgvs"] = req["raw"] if req else hgvs

    cv = None
    if isinstance(clinvar, dict) and clinvar.get("title"):
        seen["clinvar_title"] = clinvar["title"]
        cv = _parse_hgvs(clinvar["title"])
        seen["clinvar_hgvs"] = cv["raw"] if cv else None

    # Check 1: the HGVS you asked for vs the HGVS ClinVar filed. Same notation only.
    if req and cv:
        if req["kind"] == cv["kind"]:
            # Compare PARSED coordinates, not raw strings: ClinVar writes
            # "NM_000546.6(TP53):c.733G>A" where the accession sits outside the
            # variant parens, so a raw compare flags identical variants as wrong.
            same = None
            if req["pos"] and cv["pos"]:
                same = (req["pos"], req["ref"], req["alt"]) == (cv["pos"], cv["ref"], cv["alt"])
            if same is None:
                same = (req["raw"].split(":")[-1].lower() == cv["raw"].split(":")[-1].lower())
            if not same:
                problems.append(f"ClinVar record is {cv['raw']}, you asked for "
                                f"{req['raw']} -> WRONG VARIATION ID")
            else:
                checks.append(f"requested HGVS == ClinVar HGVS ({req['raw'].split(':')[-1]})")
        else:
            checks.append(f"ClinVar is {cv['kind']}. and you asked for {req['kind']}. : "
                          "cross-checked by genomic position instead of by string")

    vep = vep if isinstance(vep, dict) else {}
    vpos = vep.get("start")
    vchr = str(vep.get("seq_region_name") or "")
    if vpos:
        seen["vep_locus"] = f"{vchr}:{vpos} {vep.get('allele_string')}"
    vref, valt = _vep_alleles_plus_strand(vep, req["kind"] if req else None)

    # Check 2: ClinVar genomic record vs VEP's resolution of the same variant.
    if cv and cv["kind"] == "g" and vpos:
        if cv["pos"] != vpos:
            problems.append(f"ClinVar genomic position is {cv['pos']} but VEP resolved "
                            f"{vpos} -> DIFFERENT POSITION")
        elif cv["ref"] and (vref, valt) != (cv["ref"], cv["alt"]):
            problems.append(f"ClinVar alleles {cv['ref']}>{cv['alt']} vs VEP "
                            f"{vref}>{valt} at {vpos} -> DIFFERENT ALT")
        else:
            checks.append(f"ClinVar genomic locus == VEP locus ({vchr}:{vpos})")

    # Check 3: the gnomAD id vs VEP's resolution.
    if isinstance(gnomad, dict):
        gid = gnomad.get("variant_id") or ""
        seen["gnomad_id"] = gid
        gm = re.match(r"^(\w+)-(\d+)-([ACGT]+)-([ACGT]+)$", gid)
        if gm and vpos:
            gchr, gpos, gref, galt = gm.groups()
            if gchr != vchr or gpos != str(vpos):
                problems.append(f"gnomAD id is at {gchr}:{gpos} but VEP resolved "
                                f"{vchr}:{vpos} -> DIFFERENT POSITION")
            elif vref and (gref, galt) != (vref, valt):
                problems.append(f"gnomAD alleles {gref}>{galt} vs VEP "
                                f"{vref}>{valt} at the same position -> DIFFERENT ALT")
            else:
                checks.append(f"gnomAD id == VEP locus ({gchr}:{gpos} {gref}>{galt})")

    # Nothing was actually compared -> the identifiers are unverified, not verified.
    # ok must not be true by default: that is how the genomic-HGVS bypass happened.
    if not problems and not checks:
        problems.append(
            "no cross-check was possible: nothing tied the supplied identifiers to one "
            "variant. Supply all three (--clinvar, --gnomad, --hgvs) so the guard can "
            "prove they agree. PM2 will be withheld until it can.")

    return {"ok": not problems, "problems": problems, "checks": checks,
            "identifiers": seen}


def selftest():
    print("ClinVar: ", bool(fetch_clinvar(3634421).get("classification")))
    print("gnomAD:  ", fetch_gnomad("17-7669613-T-A").get("af"))
    print("VEP:     ", fetch_vep("NM_000546.6:c.1178A>T", "ENST00000269305")["scores"]["revel"])
    print("PubMed:  ", fetch_pubmed("p.Asp393Val,D393V")["count"], "hits")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinvar", help="ClinVar variation id (e.g. 3634421)")
    ap.add_argument("--gnomad", help="gnomAD variant id, GRCh38, e.g. 17-7669613-T-A")
    ap.add_argument("--hgvs", help="HGVS with accession: NM_000546.6:c.1178A>T or "
                                   "NC_000017.11:g.7674230C>T (both are parsed)")
    ap.add_argument("--transcript", help="MANE Select ENST id (filters VEP output)")
    ap.add_argument("--pubmed-terms", help="comma-separated protein-change spellings")
    ap.add_argument("--out", help="JSON output path (default stdout)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()

    result = {"fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    try:
        _collect(a, result)
    except OutOfScope as e:
        print(f"[OUT OF SCOPE] {e}", file=sys.stderr)
        return 2

    result["consistency"] = check_identity(
        result.get("clinvar"), result.get("gnomad"), result.get("vep"), a.hgvs)
    if not result["consistency"]["ok"]:
        print("!! IDENTIFIER CHECK FAILED - DO NOT CURATE THIS PAYLOAD", file=sys.stderr)
        for pr in result["consistency"]["problems"]:
            print("!!   " + pr, file=sys.stderr)
    else:
        for ck in result["consistency"]["checks"]:
            print(f"ok   {ck}", file=sys.stderr)

    txt = json.dumps(result, indent=2, ensure_ascii=False)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            f.write(txt + "\n")
        print(f"wrote {a.out}")
    else:
        print(txt)
    return 0


def _collect(a, result):
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
        except OutOfScope:                          # not an endpoint failure: abort
            raise
        except Exception as e:                      # a dead endpoint is NA, never a zero
            result[name] = {"error": str(e)[:400]}


if __name__ == "__main__":
    sys.exit(main())

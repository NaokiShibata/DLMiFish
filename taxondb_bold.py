from __future__ import annotations

import gzip
import hashlib
import json
import re
import time
from typing import Any, Dict, Iterable, List, Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


DEFAULT_BASE_URL = "https://portal.boldsystems.org/api"
DEFAULT_TIMEOUT_SEC = 60.0
DEFAULT_RETRIES = 3
DEFAULT_BACKOFF_SEC = 1.5
DEFAULT_USER_AGENT = "TaxonDBBuilder/0.1 (+https://github.com/NaokiShibata/TaxonDBBuilder)"
MAX_DOCUMENT_COUNT = 1_000_000


class BoldApiError(RuntimeError):
    pass


def get_bold_runtime_config(raw_cfg: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    cfg = dict(raw_cfg or {})
    base_url = str(cfg.get("base_url", DEFAULT_BASE_URL)).rstrip("/")
    timeout_sec = float(cfg.get("timeout_sec", DEFAULT_TIMEOUT_SEC))
    retries = int(cfg.get("retries", DEFAULT_RETRIES))
    backoff_sec = float(cfg.get("backoff_sec", DEFAULT_BACKOFF_SEC))
    user_agent = str(cfg.get("user_agent", DEFAULT_USER_AGENT)).strip() or DEFAULT_USER_AGENT

    if timeout_sec <= 0:
        raise BoldApiError("bold.timeout_sec must be > 0.")
    if retries < 0:
        raise BoldApiError("bold.retries must be >= 0.")
    if backoff_sec < 0:
        raise BoldApiError("bold.backoff_sec must be >= 0.")

    return {
        "base_url": base_url,
        "timeout_sec": timeout_sec,
        "retries": retries,
        "backoff_sec": backoff_sec,
        "user_agent": user_agent,
    }


def build_taxon_query(scientific_name: str) -> str:
    return f"tax:{scientific_name}"


def _parse_json_payload(text: str) -> Any:
    stripped = text.lstrip("\ufeff").strip()
    if not stripped:
        raise json.JSONDecodeError("Empty response", text, 0)

    try:
        return json.loads(stripped)
    except json.JSONDecodeError:
        pass

    decoder = json.JSONDecoder()
    values: List[Any] = []
    idx = 0
    length = len(stripped)
    while idx < length:
        while idx < length and stripped[idx].isspace():
            idx += 1
        if idx >= length:
            break
        value, next_idx = decoder.raw_decode(stripped, idx)
        values.append(value)
        idx = next_idx

    if values:
        return values if len(values) > 1 else values[0]

    line_values: List[Any] = []
    for line in stripped.splitlines():
        line = line.strip()
        if not line:
            continue
        line_values.append(json.loads(line))
    if line_values:
        return line_values if len(line_values) > 1 else line_values[0]

    raise json.JSONDecodeError("Unable to parse JSON payload", text, 0)


def _request_json(url: str, runtime_cfg: Dict[str, Any]) -> Any:
    retries = int(runtime_cfg["retries"])
    timeout_sec = float(runtime_cfg["timeout_sec"])
    backoff_sec = float(runtime_cfg["backoff_sec"])
    headers = {
        "Accept": "application/json",
        "Accept-Encoding": "gzip",
        "User-Agent": str(runtime_cfg["user_agent"]),
    }

    for attempt in range(retries + 1):
        request = Request(url, headers=headers, method="GET")
        try:
            with urlopen(request, timeout=timeout_sec) as response:
                body = response.read()
                encoding = response.headers.get("Content-Encoding", "")
                if "gzip" in encoding.lower():
                    body = gzip.decompress(body)
                charset = response.headers.get_content_charset() or "utf-8"
                text = body.decode(charset, errors="replace")
                return _parse_json_payload(text)
        except HTTPError as exc:
            if attempt >= retries:
                raise BoldApiError(f"BOLD HTTP error {exc.code}: {exc.reason}") from exc
        except URLError as exc:
            if attempt >= retries:
                raise BoldApiError(f"BOLD network error: {exc.reason}") from exc
        except json.JSONDecodeError as exc:
            raise BoldApiError(f"BOLD response JSON decode error: {exc}") from exc

        if backoff_sec > 0:
            time.sleep(backoff_sec * (attempt + 1))

    raise BoldApiError("BOLD request failed after retries.")


def _build_url(base_url: str, path: str, params: Optional[Dict[str, str]] = None) -> str:
    url = f"{base_url}{path}"
    if params:
        return f"{url}?{urlencode(params)}"
    return url


def _normalized_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def _find_values_by_key(payload: Any, normalized_targets: Iterable[str]) -> List[Any]:
    targets = set(normalized_targets)
    found: List[Any] = []

    def visit(node: Any) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if _normalized_key(str(key)) in targets:
                    found.append(value)
                visit(value)
        elif isinstance(node, list):
            for item in node:
                visit(item)

    visit(payload)
    return found


def _first_scalar(payload: Any, keys: Iterable[str]) -> Optional[str]:
    normalized_targets = [_normalized_key(key) for key in keys]
    for value in _find_values_by_key(payload, normalized_targets):
        if isinstance(value, str) and value.strip():
            return value.strip()
        if isinstance(value, (int, float)):
            return str(value)
        if isinstance(value, list):
            for item in value:
                if isinstance(item, str) and item.strip():
                    return item.strip()
                if isinstance(item, (int, float)):
                    return str(item)
    return None


def _extract_successful_terms(preprocess_payload: Any) -> List[str]:
    if isinstance(preprocess_payload, dict):
        successful_terms = preprocess_payload.get("successful_terms")
    else:
        successful_terms = None

    if not isinstance(successful_terms, list):
        return []

    terms: List[str] = []
    for term in successful_terms:
        if isinstance(term, str) and term.strip():
            terms.append(term.strip())
            continue
        if not isinstance(term, dict):
            continue
        matched = term.get("matched")
        if isinstance(matched, str) and matched.strip():
            terms.append(matched.strip())
            continue
        scope = term.get("scope")
        field = term.get("field")
        value = term.get("value")
        if isinstance(scope, str) and isinstance(value, str):
            if isinstance(field, str) and field.strip():
                terms.append(f"{scope}:{field}:{value}")
            else:
                terms.append(f"{scope}:{value}")
    return [term for term in terms if term]


def preprocess_query(query: str, runtime_cfg: Dict[str, Any]) -> Tuple[str, Any]:
    url = _build_url(runtime_cfg["base_url"], "/query/preprocessor", {"query": query})
    payload = _request_json(url, runtime_cfg)
    successful_terms = _extract_successful_terms(payload)
    normalized_query = ";".join(successful_terms) if successful_terms else query
    return normalized_query, payload


def fetch_specimen_count(normalized_query: str, runtime_cfg: Dict[str, Any]) -> Optional[int]:
    url = _build_url(
        runtime_cfg["base_url"],
        "/summary",
        {
            "query": normalized_query,
            "fields": "specimens",
            "reduce_operation": "count",
        },
    )
    payload = _request_json(url, runtime_cfg)
    count_value = None
    if isinstance(payload, dict):
        counts = payload.get("counts")
        if isinstance(counts, dict):
            count_value = counts.get("specimens")
        if count_value is None:
            count_value = payload.get("specimens")
        if count_value is None:
            count_value = payload.get("count")
    if count_value is None:
        return None
    if isinstance(count_value, list) and count_value:
        count_value = count_value[0]
    try:
        return int(count_value)
    except (TypeError, ValueError):
        return None


def submit_query(normalized_query: str, runtime_cfg: Dict[str, Any], extent: str = "full") -> Tuple[str, Any]:
    url = _build_url(runtime_cfg["base_url"], "/query", {"query": normalized_query, "extent": extent})
    payload = _request_json(url, runtime_cfg)
    query_id = None
    if isinstance(payload, dict):
        for key in ("query_id", "queryId", "id", "token"):
            value = payload.get(key)
            if isinstance(value, str) and value.strip():
                query_id = value.strip()
                break
    if not query_id:
        raise BoldApiError("BOLD query response did not contain a query_id.")
    return query_id, payload


def download_documents(query_id: str, runtime_cfg: Dict[str, Any], fmt: str = "json") -> Any:
    encoded_query_id = quote(query_id, safe="")
    url = _build_url(runtime_cfg["base_url"], f"/documents/{encoded_query_id}/download", {"format": fmt})
    return _request_json(url, runtime_cfg)


def extract_document_rows(documents_payload: Any) -> List[Dict[str, Any]]:
    if isinstance(documents_payload, list):
        return [row for row in documents_payload if isinstance(row, dict)]
    if isinstance(documents_payload, dict):
        for key in ("documents", "records", "items", "data", "results"):
            value = documents_payload.get(key)
            if isinstance(value, list):
                return [row for row in value if isinstance(row, dict)]
        if any(_first_scalar(documents_payload, [key]) for key in ("processid", "sampleid", "marker_code", "nucleotides")):
            return [documents_payload]
    return []


def parse_accession_tokens(raw_value: Optional[str]) -> List[str]:
    if not raw_value:
        return []
    tokens = re.split(r"[,\s;|/]+", raw_value.strip())
    cleaned: List[str] = []
    for token in tokens:
        value = token.strip()
        if value and value not in cleaned:
            cleaned.append(value)
    return cleaned


def _normalize_marker_text(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def _marker_matches(marker_text: str, candidate: str, exact_only: bool) -> bool:
    left = _normalize_marker_text(marker_text)
    right = _normalize_marker_text(candidate)
    if not left or not right:
        return False
    if left == right:
        return True
    if exact_only:
        return False
    return left.startswith(right) or right.startswith(left) or right in left or left in right


def resolve_marker_from_record(
    marker_text: Optional[str],
    marker_keys: List[str],
    marker_map: Dict[str, Dict[str, Any]],
) -> Optional[Tuple[str, str]]:
    if not marker_text:
        return None

    strategies = ("marker_codes", "aliases", "phrases", "key")
    for strategy in strategies:
        for marker_key in marker_keys:
            marker_cfg = marker_map.get(marker_key, {})
            if strategy == "marker_codes":
                candidates = list((marker_cfg.get("bold") or {}).get("marker_codes") or [])
                exact_only = True
            elif strategy == "aliases":
                candidates = list(marker_cfg.get("aliases") or [])
                exact_only = False
            elif strategy == "phrases":
                candidates = list(marker_cfg.get("phrases") or [])
                exact_only = False
            else:
                candidates = [marker_key]
                exact_only = False
            for candidate in candidates:
                if not isinstance(candidate, str) or not candidate.strip():
                    continue
                if _marker_matches(marker_text, candidate, exact_only=exact_only):
                    return marker_key, marker_text
    return None


def _sanitize_sequence(sequence: Optional[str]) -> str:
    if not sequence:
        return ""
    cleaned = re.sub(r"[^A-Za-z]", "", sequence).upper()
    return cleaned


def normalize_bold_row(
    raw_row: Dict[str, Any],
    marker_keys: List[str],
    marker_map: Dict[str, Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    marker_text = _first_scalar(raw_row, ["marker_code", "marker", "markercode"])
    marker_match = resolve_marker_from_record(marker_text, marker_keys, marker_map)
    if not marker_match:
        return None

    sequence = _sanitize_sequence(
        _first_scalar(
            raw_row,
            [
                "nucleotides",
                "nuc",
                "sequence",
                "nucleotide",
                "nucleotidesequence",
            ],
        )
    )
    if not sequence:
        return None

    processid = _first_scalar(raw_row, ["processid", "process_id"])
    sampleid = _first_scalar(raw_row, ["sampleid", "sample_id"])
    accession_raw = _first_scalar(raw_row, ["insdcacs", "insdc_acs", "genbank_accession"])
    taxon_name = _first_scalar(
        raw_row,
        [
            "species",
            "scientific_name",
            "taxon_name",
            "identification",
            "taxonomy_species",
        ],
    )

    marker_key, marker_label = marker_match
    source_record_id = processid or sampleid or accession_raw
    if not source_record_id:
        seed = f"{taxon_name or ''}|{marker_key}|{sequence}"
        source_record_id = f"bold_{hashlib.sha1(seed.encode('utf-8')).hexdigest()[:16]}"

    return {
        "source_record_id": source_record_id,
        "processid": processid,
        "sampleid": sampleid,
        "accession": accession_raw,
        "taxon_name": taxon_name,
        "marker_key": marker_key,
        "marker_label": marker_label,
        "sequence": sequence,
        "raw_row": raw_row,
    }


def fetch_bold_records_for_taxon(
    scientific_name: str,
    marker_keys: List[str],
    marker_map: Dict[str, Dict[str, Any]],
    bold_cfg: Optional[Dict[str, Any]] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    runtime_cfg = get_bold_runtime_config(bold_cfg)
    raw_query = build_taxon_query(scientific_name)
    normalized_query, preprocess_payload = preprocess_query(raw_query, runtime_cfg)
    specimen_count = fetch_specimen_count(normalized_query, runtime_cfg)
    if specimen_count is not None and specimen_count > MAX_DOCUMENT_COUNT:
        raise BoldApiError(
            f"BOLD query exceeds maximum downloadable records ({specimen_count} > {MAX_DOCUMENT_COUNT})."
        )
    if specimen_count == 0:
        return [], {
            "raw_query": raw_query,
            "normalized_query": normalized_query,
            "specimen_count": 0,
            "preprocess_payload": preprocess_payload,
            "query_id": None,
            "downloaded_rows": 0,
            "matched_rows": 0,
        }

    query_id, query_payload = submit_query(normalized_query, runtime_cfg, extent="full")
    documents_payload = download_documents(query_id, runtime_cfg, fmt="json")
    rows = extract_document_rows(documents_payload)

    normalized_rows: List[Dict[str, Any]] = []
    for row in rows:
        normalized = normalize_bold_row(row, marker_keys, marker_map)
        if normalized:
            normalized_rows.append(normalized)

    stats = {
        "raw_query": raw_query,
        "normalized_query": normalized_query,
        "specimen_count": specimen_count,
        "preprocess_payload": preprocess_payload,
        "query_payload": query_payload,
        "query_id": query_id,
        "downloaded_rows": len(rows),
        "matched_rows": len(normalized_rows),
    }
    return normalized_rows, stats

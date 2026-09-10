"""
registration.py — is the constraint file actually pre-registered, and does it still match?

THE GATE THIS IMPLEMENTS
------------------------
Anatomy_Test.md, step 2:

    Write and timestamp the constraint file ... Commit it, hash it, and register it
    publicly — OSF takes minutes and is free.
    **Gate: registration timestamp precedes every result file.**

and step 6's risk table names circularity first:

    If a constraint was learned from deconvolution output, the test is rigged and a
    reviewer will find it. → Every constraint cites histology or in-situ hybridization,
    and the file is timestamped before any output exists.

Nothing checked either. The constraint file was hashed, which proves it has not changed
since *some* moment, but not that the moment came before the results. A hash with no
timestamp is a checksum, not a pre-registration, and the difference is the whole
credibility of the design: "we predicted this" versus "we can show we predicted this".

WHAT THIS CAN AND CANNOT ESTABLISH
----------------------------------
It cannot create trust that does not exist. A local file claiming a registration date is
worth nothing to a reviewer — anyone can write one. What it *can* do is:

  * refuse to let the project describe itself as pre-registered when it is not;
  * detect the case that would actually be fatal — the constraint file changing after
    registration — by comparing the recorded hash against the live one;
  * check the ordering the protocol asks for, that results were produced after the
    registration and not before.

So the honest default is `UNREGISTERED`, stated in the report and in RESULTS.md, until
someone registers the hash publicly and records the receipt here.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path

from ivygap import config
from ivygap.anatomic import constraints as K

#: Where the public-registration receipt lives once one exists. Committed, so the claim
#: travels with the repository rather than with one machine.
REGISTRATION_PATH = config.PROJECT_ROOT / "REGISTRATION.json"

#: The fields a receipt must carry to mean anything.
REQUIRED_FIELDS = ("registry", "url", "registered_utc", "constraint_freeze_hash")


#: HTTP meanings for an OSF registration URL, learned the hard way.
#:
#: On 2026-09-10 this project recorded a receipt for https://osf.io/nj78m, reported
#: REGISTERED, and embedded that in every artefact a run produced. The registration had
#: been deleted when its approval was rejected — the URL returned 410 — and nothing
#: noticed, because the checker only ever compared two hashes and never asked whether the
#: thing it pointed at still existed. A receipt for a registration that is gone is worse
#: than no receipt: it is a false claim of pre-registration, made automatically, in
#: writing, in every result file.
LIVENESS = {
    200: ("REGISTERED", "public and readable"),
    401: ("PENDING", "exists but is not public yet — submitted and awaiting approval, "
                     "or embargoed. The submission timestamp is already fixed; what is "
                     "missing is public readability, which is the point of registering."),
    403: ("PENDING", "exists but access is restricted"),
    404: ("UNREACHABLE", "no such registration at that URL"),
    410: ("WITHDRAWN", "the registration was withdrawn or deleted. A recorded receipt "
                       "pointing at it is a false claim and must not be reported as "
                       "REGISTERED."),
}


def check_url_live(url: str, timeout: float = 15.0) -> tuple[str | None, str]:
    """
    Ask OSF whether the recorded registration actually exists and is public.

    Returns (state_or_None, explanation). `None` means the question could not be
    answered — offline, DNS down, a non-OSF registry — and an unanswerable question must
    never silently downgrade a real registration, so callers keep the recorded state and
    say the check did not run.
    """
    import re
    import ssl
    import subprocess
    import urllib.error
    import urllib.request

    m = re.search(r"osf\.io/([a-z0-9]{5})", url or "", re.I)
    if not m:
        return None, f"not an OSF registration URL, so liveness was not checked: {url!r}"

    api = f"https://api.osf.io/v2/registrations/{m.group(1)}/"

    # Verify certificates properly. This machine's Python has no usable default trust
    # store, so `urlopen` raises URLError on every https call — which would make the
    # liveness check silently unable to run, which is the failure it exists to prevent.
    ctx = None
    try:
        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
    except Exception:
        try:
            ctx = ssl.create_default_context()
        except Exception:
            ctx = None

    def _fetch() -> int | None:
        try:
            req = urllib.request.Request(
                api, headers={"Accept": "application/vnd.api+json"})
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
                return r.status
        except urllib.error.HTTPError as e:
            return e.code
        except Exception:
            return None

    code = _fetch()
    if code is None:
        # Fall back to curl, which carries the system trust store on macOS.
        try:
            out = subprocess.run(
                ["curl", "-s", "-o", "/dev/null", "-w", "%{http_code}",
                 "--max-time", str(int(timeout)),
                 "-H", "Accept: application/vnd.api+json", api],
                capture_output=True, text=True, timeout=timeout + 5)
            code = int(out.stdout.strip()) or None
        except Exception:
            code = None
    if code is None:
        return None, "liveness check could not run (no network); recorded state kept"

    state, why = LIVENESS.get(code, (None, f"unexpected HTTP {code}"))
    return state, f"HTTP {code} — {why}"


@dataclass
class RegistrationStatus:
    state: str                     # REGISTERED | PENDING | UNREGISTERED |
                                   # HASH_MISMATCH | MALFORMED | UNREACHABLE | WITHDRAWN
    freeze_hash: str
    recorded_hash: str | None = None
    registered_utc: str | None = None
    registry: str | None = None
    url: str | None = None
    results_predating_registration: list = None
    verdict: str = ""

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)

    @property
    def is_registered(self) -> bool:
        return self.state == "REGISTERED"


def _result_files() -> list[Path]:
    out = []
    for d in (config.ANATOMIC_DIR, config.BENCH_DIR, config.SURVIVAL_DIR):
        if d.exists():
            out.extend(p for p in d.rglob("*") if p.is_file())
    return out


def status(check_results: bool = True, check_live: bool = True) -> RegistrationStatus:
    """
    Read the registration receipt and check it against the live constraint file.

    `check_results=False` skips the file-ordering check, which is what a run wants
    *before* it has written anything. `check_live=False` skips the network call that
    asks OSF whether the recorded registration still exists.
    """
    live = K.freeze_hash()

    if not REGISTRATION_PATH.exists():
        return RegistrationStatus(
            state="UNREGISTERED", freeze_hash=live,
            verdict=(
                "UNREGISTERED: the constraint file is hashed and committed, but no "
                "public registration receipt exists. The hash proves the constraints "
                "have not changed; it does not prove they were written before the "
                "results, which is what the pre-registration claim rests on. Register "
                f"{live[:16]}... publicly (OSF takes minutes) and record the receipt in "
                f"{REGISTRATION_PATH.name}. Until then this project must not describe "
                "itself as pre-registered."),
        )

    try:
        rec = json.loads(REGISTRATION_PATH.read_text())
    except Exception as exc:                                # noqa: BLE001
        return RegistrationStatus(
            state="MALFORMED", freeze_hash=live,
            verdict=f"MALFORMED: {REGISTRATION_PATH.name} is not readable JSON ({exc}).")

    missing = [f for f in REQUIRED_FIELDS if not rec.get(f)]
    if missing:
        return RegistrationStatus(
            state="MALFORMED", freeze_hash=live,
            verdict=f"MALFORMED: {REGISTRATION_PATH.name} is missing {missing}.")

    recorded = str(rec["constraint_freeze_hash"])
    if recorded != live:
        return RegistrationStatus(
            state="HASH_MISMATCH", freeze_hash=live, recorded_hash=recorded,
            registered_utc=rec.get("registered_utc"), registry=rec.get("registry"),
            url=rec.get("url"),
            verdict=(
                f"HASH MISMATCH — this is the fatal one. The registered constraint file "
                f"hashes to {recorded[:16]}... and the live one to {live[:16]}.... The "
                f"constraints have been edited since registration, so every result "
                f"computed under the live file is outside the pre-registration. Either "
                f"restore the registered constraints or register the new ones and treat "
                f"everything before as superseded."),
        )

    # The recorded hash matches. Before claiming REGISTERED, ask whether the thing the
    # receipt points at still exists — a withdrawn registration must never be reported
    # as a live one. See LIVENESS: this project asserted REGISTERED against a deleted
    # URL for a full day because nothing ever asked.
    if check_live:
        live_state, why = check_url_live(str(rec.get("url") or ""))
        if live_state in ("WITHDRAWN", "UNREACHABLE"):
            return RegistrationStatus(
                state=live_state, freeze_hash=live, recorded_hash=recorded,
                registered_utc=rec.get("registered_utc"), registry=rec.get("registry"),
                url=rec.get("url"),
                verdict=(
                    f"{live_state}: {REGISTRATION_PATH.name} records "
                    f"{rec.get('url')}, but that registration {why}. The receipt is "
                    f"stale and the project is NOT registered. Do not report these "
                    f"results as pre-registered until a live registration is recorded."),
            )
        if live_state == "PENDING":
            return RegistrationStatus(
                state="PENDING", freeze_hash=live, recorded_hash=recorded,
                registered_utc=rec.get("registered_utc"), registry=rec.get("registry"),
                url=rec.get("url"),
                verdict=(
                    f"PENDING: {rec.get('url')} exists and its submission timestamp is "
                    f"fixed, but it is not publicly readable ({why}). The constraints "
                    f"are committed and unchanged. Until it is public a reader cannot "
                    f"verify the claim, so results may be described as SUBMITTED for "
                    f"registration, not as pre-registered."),
            )


    predating = []
    if check_results:
        try:
            reg_at = datetime.fromisoformat(
                str(rec["registered_utc"]).replace("Z", "+00:00"))
            if reg_at.tzinfo is None:
                reg_at = reg_at.replace(tzinfo=timezone.utc)
            for p in _result_files():
                mtime = datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc)
                if mtime < reg_at:
                    predating.append(str(p.relative_to(config.PROJECT_ROOT)))
        except Exception:                                   # noqa: BLE001
            pass

    if predating:
        return RegistrationStatus(
            state="REGISTERED", freeze_hash=live, recorded_hash=recorded,
            registered_utc=rec.get("registered_utc"), registry=rec.get("registry"),
            url=rec.get("url"), results_predating_registration=sorted(predating),
            verdict=(
                f"REGISTERED, but {len(predating)} result file(s) are OLDER than the "
                f"registration. The protocol's gate is that the registration precedes "
                f"every result. Those files were produced before the constraints were "
                f"registered and must be regenerated or labelled superseded."),
        )

    return RegistrationStatus(
        state="REGISTERED", freeze_hash=live, recorded_hash=recorded,
        registered_utc=rec.get("registered_utc"), registry=rec.get("registry"),
        url=rec.get("url"), results_predating_registration=[],
        verdict=(f"REGISTERED: constraint file {live[:16]}... registered at "
                 f"{rec.get('registered_utc')} with {rec.get('registry')}, and every "
                 f"result file postdates it."),
    )


def write_template(path: Path | None = None) -> Path:
    """
    Write a receipt template with the CURRENT hash filled in, for someone to complete
    after registering publicly. Never marks the project registered on its own.
    """
    path = path or (config.PROJECT_ROOT / "REGISTRATION.template.json")
    path.write_text(json.dumps({
        "_instructions": (
            "Register the constraint_freeze_hash below on a public timestamping "
            "registry (OSF, Zenodo, or an OpenTimestamps proof), then fill in registry, "
            "url and registered_utc and rename this file to REGISTRATION.json. Do not "
            "fill it in without registering: the point is the external record, and a "
            "self-asserted date is worth nothing to a reviewer."),
        "registry": "",
        "url": "",
        "registered_utc": "",
        "constraint_freeze_hash": K.freeze_hash(),
    }, indent=2))
    return path

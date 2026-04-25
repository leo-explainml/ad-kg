"""FDA FAERS ingest: ROR computation for AD protective signals."""

from __future__ import annotations

import logging
import math
from typing import Any

import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from ad_kg.models import FAERSReport

logger = logging.getLogger(__name__)

_OPENFDA_BASE = "https://api.fda.gov/drug/event.json"

_DEFAULT_DRUGS = [
    # GLP-1 agonists
    "semaglutide",
    "liraglutide",
    "exenatide",
    "dulaglutide",
    "tirzepatide",
    "albiglutide",    # appeared in gwas_snp_to_drug via GLP1R
    "lixisenatide",   # appeared in gwas_snp_to_drug via GLP1R
    # Biguanides
    "metformin",
    # SGLT2 inhibitors
    "empagliflozin",
    "dapagliflozin",
    "canagliflozin",
    "ertugliflozin",
    # Thiazolidinediones — PPARG agonists with GWAS bridge gene support and AD biology
    "pioglitazone",
]

# MedDRA terms associated with Alzheimer's / cognitive impairment
_AD_REACTIONS = [
    "Dementia",
    "Alzheimer's disease",
    "Memory impairment",
    "Cognitive disorder",
    "Dementia Alzheimer's type",
    "Cognitive impairment",
]

# Cohorts for sensitivity analysis.
# Each value is an OpenFDA Lucene filter appended to the base search query.
# "all" = no filter (baseline ROR); others restrict the reporting population.
_COHORTS: dict[str, str | None] = {
    "all": None,
    "t2dm": 'patient.drug.drugindication:"type 2 diabetes mellitus"',
    "elderly": "patient.patientagegroup:6",
    "post_2020": "receivedate:[20200101 TO 20991231]",
}


@retry(
    retry=retry_if_exception_type((requests.RequestException, IOError)),
    wait=wait_exponential(multiplier=1, min=2, max=60),
    stop=stop_after_attempt(5),
    reraise=True,
)
def _fda_get(params: dict[str, Any]) -> dict[str, Any]:
    """GET the OpenFDA drug/event endpoint with retry."""
    resp = requests.get(_OPENFDA_BASE, params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()


def compute_ror(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    """Compute Reporting Odds Ratio with 95% CI using the log method.

    2x2 contingency table:
        a = drug+AD reports
        b = drug+notAD reports
        c = notdrug+AD reports
        d = notdrug+notAD reports

    Returns:
        (ror, ci_lower, ci_upper)
    """
    # Add 0.5 continuity correction for zeros
    a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    log_ror = math.log(a_ * d_) - math.log(b_ * c_)
    se = math.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    z = 1.96
    log_lower = log_ror - z * se
    log_upper = log_ror + z * se
    return math.exp(log_ror), math.exp(log_lower), math.exp(log_upper)


def _count_reports(
    drug: str,
    reaction: str | None = None,
    cohort_filter: str | None = None,
) -> int:
    """Count FDA adverse event reports for a drug (+ optional reaction and cohort filter)."""
    search_parts = [f'patient.drug.medicinalproduct:"{drug}"']
    if reaction:
        search_parts.append(f'patient.reaction.reactionmeddrapt:"{reaction}"')
    if cohort_filter:
        search_parts.append(cohort_filter)
    params: dict[str, Any] = {
        "search": " AND ".join(search_parts),
        "limit": 1,
    }
    try:
        data = _fda_get(params)
        return data.get("meta", {}).get("results", {}).get("total", 0)
    except Exception as exc:
        logger.debug("FDA count failed for drug=%r reaction=%r: %s", drug, reaction, exc)
        return 0


def _count_all_drug_reports(drug: str, cohort_filter: str | None = None) -> int:
    """Total reports for a drug regardless of reaction."""
    return _count_reports(drug, cohort_filter=cohort_filter)


def _count_total_reports_with_reaction(
    reaction: str, cohort_filter: str | None = None
) -> int:
    """Total reports with a given reaction regardless of drug."""
    search_parts = [f'patient.reaction.reactionmeddrapt:"{reaction}"']
    if cohort_filter:
        search_parts.append(cohort_filter)
    params: dict[str, Any] = {
        "search": " AND ".join(search_parts),
        "limit": 1,
    }
    try:
        data = _fda_get(params)
        return data.get("meta", {}).get("results", {}).get("total", 0)
    except Exception as exc:
        logger.debug("FDA reaction count failed for %r: %s", reaction, exc)
        return 0


def _count_total_reports(cohort_filter: str | None = None) -> int:
    """Approximate total number of reports in the database (optionally filtered to a cohort)."""
    params: dict[str, Any] = {"limit": 1}
    if cohort_filter:
        params["search"] = cohort_filter
    try:
        data = _fda_get(params)
        return data.get("meta", {}).get("results", {}).get("total", 10_000_000)
    except Exception:
        return 10_000_000


def ingest_faers(
    drug_list: list[str] | None = None,
    cohorts: dict[str, str | None] | None = None,
) -> list[FAERSReport]:
    """Ingest FAERS data and compute ROR for AD-related reactions across cohorts.

    Args:
        drug_list: Drug names to query. Defaults to GLP-1/metformin/SGLT2 list.
        cohorts: Dict of cohort_name -> OpenFDA search filter. Defaults to
            _COHORTS (all, t2dm, elderly, post_2020).

    Returns:
        List of FAERSReport objects with computed ROR, one per drug/reaction/cohort.
    """
    if drug_list is None:
        drug_list = _DEFAULT_DRUGS
    if cohorts is None:
        cohorts = _COHORTS

    logger.info("Fetching FAERS for %d drugs × %d cohorts", len(drug_list), len(cohorts))
    reports: list[FAERSReport] = []

    for cohort_name, cohort_filter in cohorts.items():
        logger.info("Cohort: %s", cohort_name)
        total_reports = _count_total_reports(cohort_filter)

        for drug in drug_list:
            logger.info("  drug=%r cohort=%s", drug, cohort_name)
            drug_total = _count_all_drug_reports(drug, cohort_filter)

            for reaction in _AD_REACTIONS:
                a = _count_reports(drug, reaction, cohort_filter)
                if a == 0:
                    continue

                c = _count_total_reports_with_reaction(reaction, cohort_filter)
                b = drug_total - a
                d = total_reports - drug_total - c + a

                b = max(b, 0)
                d = max(d, 0)

                ror, ci_lower, ci_upper = compute_ror(a, b, c, d)

                reports.append(
                    FAERSReport(
                        drug_name=drug,
                        reaction=reaction,
                        ror=ror,
                        ci_lower=ci_lower,
                        ci_upper=ci_upper,
                        report_count=a,
                        cohort=cohort_name,
                    )
                )
                logger.debug(
                    "    %s / %s [%s]: a=%d ROR=%.3f [%.3f, %.3f]",
                    drug, reaction, cohort_name, a, ror, ci_lower, ci_upper,
                )

    logger.info("FAERS reports generated: %d", len(reports))
    return reports

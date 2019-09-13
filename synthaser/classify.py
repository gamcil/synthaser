#!/usr/bin/env python3

"""This module stores all of the routines for classification of Synthase objects into
biosynthetic categories.
"""


PKS_TYPES = {
    "Type I PKS": ["PKS", "PKS_KS"],
    "Type II PKS": ["CLF", "KAS_I_II"],
    "Type III PKS": ["CHS_like", "KAS_III"],
}


def assign_PKS_type(domains):
    """Classify a PKS as Type I, II or III based on their KS CDD domain name.

    If the supplied domains contain 2 or more ketoacyl-synthase (KS) domains, the PKS
    will be classified as multi-modular.

    Otherwise, it will be classified based on the following rules.

    Type I:
        KS or PKS_KS
    Type II:
        CLF or KAS_I_II
    Type III:
        CHS_like

    Note that all of these domains are stored as KS domains. During filtering, those
    with the maximum score are chosen; Type I PKSs generally will have a higher score
    for KS/PKS_KS CDs, type II PKSs for CLF/KAS_I_II, and type III PKSs for CHS_like
    (i.e. chalcone synthase like).

    This function references the rules stored in `classify.PKS_TYPES`. Thus, if you want
    to add additional types or different conserved domain names, simply update that
    dictionary. Any changes will be reflected in the output of this function.

    Parameters
    ----------
    domains : list
        Domain objects in a Synthase.

    Returns
    -------
    type : str
        The PKS type of this Synthase based on given Domains (see above).

    Raises
    ------
    ValueError
        If no type could be assigned.
    """
    domain_types = [domain.type for domain in domains]
    if domain_types.count("KS") >= 2:
        return "multi-modular PKS"
    ks = domains[domain_types.index("KS")].domain
    for type, conserved_domains in PKS_TYPES.items():
        if ks in conserved_domains:
            return type
    raise ValueError("Failed to assign PKS type")


def assign_T1PKS_subtype(domains):
    """Classify a Type I PKS as highly (HR), partially (PR) or non-reducing (NR), or other
    (PKS-like).

    These classifications are based on the presence or absence of reducing domains, i.e.
    enoyl-reductases (ER), keto-reductases (KR) and dehydratases (DH).

    HR-PKS:
        ER, KR and DH
    PR-PKS:
        Any, but not all, of the reduction domains.
    NR-PKS:
        At least KS and acyltransferase (AT) domains; if the function reaches this
        point, any synthase with reducing domains should have already returned. Thus,
        we just need to check that we still have a full PKS module.
    PKS-like:
        The remainder; usually just anything with a KS domain.

    Parameters
    ----------
    domains : list
        List of domain types in a Synthase.

    Returns
    -------
    str
        Sub-classification based on given domain types (see above).
    """
    if {"ER", "KR", "DH"}.issubset(domains):
        return "HR-PKS"
    if any(domain in domains for domain in {"ER", "KR", "DH"}):
        return "PR-PKS"
    if {"KS", "AT"}.issubset(domains):
        return "NR-PKS"
    return "PKS-like"


def assign_NRPS_type(domains):
    """Classify an NRPS as full (NRPS) or partial (NRPS-like).

    NRPS:
        At least one full module, containing adenylation (A), thiolation (T) and
        condensation (C) domains.
    NRPS-like:
        The remainder; anything with an A domain.

    Parameters
    ----------
    domains : list
        List of domain types in a Synthase.

    Returns
    -------
    str
        NRPS subclassification based on given domain types (see above).
    """
    if {"A", "ACP", "C"}.issubset(domains) or {"A", "T", "C"}.issubset(domains):
        return "NRPS"
    return "NRPS-like"


def assign_broad_type(domains):
    """Classify a Synthase as PKS, NRPS or Hybrid PKS-NRPS based on its Domains.

    Hybrid (PKS-NRPS):
        Both beta-ketoacyl synthase (KS) and adenylation (A) domains
    Polyketide synthase (PKS):
        KS domain
    Nonribosomal peptide synthase (NRPS):
        A domain

    Parameters
    ----------
    synthase : Synthase
        A Synthase object with domain hits.

    Returns
    -------
    str
        The biosynthetic type of the given Synthase (hybrid, pks, nrps).

    Raises
    ------
    ValueError
        If no identifying domain (KS, A) is found.
    """
    types = [domain.type for domain in domains]
    if {"KS", "A"}.issubset(types) or {"KS", "C"}.issubset(types):
        return "Hybrid"
    if "KS" in types:
        return assign_PKS_type(domains)
    if "A" in types:
        return assign_NRPS_type(types)
    raise ValueError("Could not find an identifying domain")


def classify_synthase(synthase):
    """Classify a Synthase.

    First, assign a broad biosynthetic type to the Synthase (PKS, NRPS, Hybrid).
    Then, further classification is performed on Type I PKSs.
    """
    synthase.type = assign_broad_type(synthase.domains)
    if synthase.type == "Type I PKS":
        synthase.subtype = assign_T1PKS_subtype(synthase.domain_types)
    else:
        synthase.subtype = synthase.type

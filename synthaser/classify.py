#!/usr/bin/env python3

"""
This module stores all of the routines for classification of Synthase objects into
biosynthetic categories.

The main function provided by this module is `classify`, which accepts a Synthase object
and runs through all the necessary classification steps for a given synthase type.

For example, if we instantiate a new `Synthase` object for a typical HR-PKS:

>>> from synthaser.models import Domain, Synthase
>>> synthase = Synthase(
...     header='sequence',
...     sequence='MKDESTMSD.....',
...     domains=[
...         Domain(type='KS', domain='PKS_KS'),
...         Domain(type='KR'),
...         Domain(type='ER'),
...         Domain(type='DH'),
...         ...
...     ]
... )

and call `classify`:

>>> from synthaser import classify
>>> classify.classify(synthase)

The function will call a series of functions to classify it.
First, it calls `assign_broad_type`, which assigns the broader biosynthetic
category ('PKS') given the first `models.Domain` instance with `type` 'KS'.
Then, it calls `assign_PKS_type`, which assigns a PKS type based on the specific
conserved domain name (`domain` attribute) of the KS domain object. These relationships
are stored (and can be altered) in the `PKS_TYPES` dictionary. In the example case,
the KS object is a 'PKS_KS' hit, and so the synthase is designated a 'Type I PKS'.

>>> synthase.type
'Type I PKS'

Finally, the synthase is further classified into a T1PKS subtype via
`assign_T1PKS_subtype`, which looks for the presence of reductive (ER, KR, DH) domains.
Since we have the full set, this synthase is a highly-reducing PKS, and as such is
assigned the subtype 'HR-PKS'.

>>> synthase.subtype
'HR-PKS'

When classifying an NRPS, the only options are NRPS (a full A, T, C module is present)
and NRPS-like.

PKS and NRPS share some related domains, which in turn will return the same set of
conserved domains when searched against the CDD. For example, an ACP domain from a PKS
and a thiolation (T) domain from an NRPS will both return 'PKS_PP', 'AcpP' and
'PP-binding' conserved domain hits. To account for this, the function
`rename_NRPS_domains` will iterate over the domain set of an NRPS or PKS-NRPS hybrid
synthase object, and rename domains based on convention.
For example, if we have an NRPS object prior to classification:

>>> synthase.architecture
A-ACP-C-A-ACP-C-TR

we can call `rename_NRPS_domains`:

>>> classify.rename_NRPS_domains(synthase)
>>> synthase.architecture
A-T-C-A-T-C-R

In this case, ACP domains are renamed T, and the thioreductase (TR) is renamed R.
"""


PKS_TYPES = {
    "SCP-x thiolase": ["SCP-x_thiolase"],
    "HMG-CoA synthase": ["HMG-CoA-S_euk"],
    "3-ketoacyl-CoA thiolase": ["thiolase", "PLN02287"],
    "FAS I/II": ["PRK07314", "KAS_I_II", "FabF", "FabB"],
    "Type I PKS": ["PKS", "PKS_KS"],
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
    ks = domains[domain_types.index("KS")].domain
    if domain_types.count("KS") >= 2:
        return "multi-modular PKS"
    for type, conserved_domains in PKS_TYPES.items():
        if ks in conserved_domains:
            return type
    if "AT" not in domain_types:
        return "trans-AT PKS"
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


# TODO: adjust for architectures with weird orders, e.g. PKS-NRPS that starts with A
#       before the PKS module ... try to specifically recognize an NRPS module
def rename_NRPS_domains(synthase):
    """Replace domain types in Hybrid and NRPS Synthases.

    The acyl carrier protein (ACP) domain in PKSs is homologous to the thioester
    domain of the peptide carrier protein (PCP) domain in NRPSs, and as such, both
    PKS and NRPS will report the same conserved domain hit. In NRPS, it is
    convention to name these T, i.e.::

        A-ACP-C --> A-T-C

    In hybrid PKS-NRPS, this replacement is made in the NRPS module of the synthase.
    Thus, this function looks for a condensation (C) domain that typically signals
    the beginning of such a module, and replaces any ACP with T after that domain.

    An example PKS-NRPS domain architecture may resemble::

        KS-AT-DH-ER-KR-ACP-C-A-T-R

    Thioester reductase (TR) domains are generally written as R in NRPS, thus the
    replacement here.

    Parameters
    ----------
    synthase : models.Synthase
        An NRPS or Hybrid Synthase object.

    Raises
    ------
    ValueError
        If `type` is not 'Hybrid' or 'NRPS'.
    """
    if synthase.type not in ("Hybrid", "NRPS"):
        raise ValueError("Expected 'Hybrid' or 'NRPS'")
    start, replace = 0, {"ACP": "T", "TR": "R"}
    if synthase.type == "Hybrid":
        for start, domain in enumerate(synthase.domains):
            if domain.type in ("C", "E"):
                break
    for domain in synthase.domains[start:]:
        if domain.type in replace:
            domain.type = replace[domain.type]


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


def classify(synthase):
    """Classify a Synthase.

    First, assign a broad biosynthetic type to the Synthase (PKS, NRPS, Hybrid).
    Then, further classification is performed on Type I PKSs.
    """
    synthase.type = assign_broad_type(synthase.domains)
    if synthase.type == "Type I PKS":
        synthase.subtype = assign_T1PKS_subtype(synthase.domain_types)
    else:
        if synthase.type in ("Hybrid", "NRPS"):
            rename_NRPS_domains(synthase)
        synthase.subtype = synthase.type

#!/usr/bin/env python3

"""
This module stores all of the routines for classification of Synthase objects into
iosynthetic categories.

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

import json
import logging


LOG = logging.getLogger(__name__)


PKS_TYPES = {
    "Thiolase": ["SCP-x_thiolase", "thiolase", "PLN02287"],
    "HMG-CoA synthase": ["HMG-CoA-S_euk"],
    "FAS": ["PRK07314", "KAS_I_II", "KAS_III", "FabF", "FabB"],
    "Type I PKS": ["PKS", "PKS_KS"],
    "Type III PKS": ["CHS_like"],
}


def tester(path):
    import json
    from synthaser.models import Domain

    with open(path) as fp:
        d = json.load(fp)
        rg = RuleGraph.from_dict(d)

    pks = [
        Domain(type="KS", domain="PKS_KS"),
        # Domain(type="KS", domain="PKS_KS_BLAH"),
        Domain(type="AT", domain="PKS_AT"),
        Domain(type="ER", domain="PKS_ER"),
    ]

    nrps = [
        Domain(type="A"),
        Domain(type="ACP"),
        Domain(type="C"),
    ]

    return rg, nrps


def _traverse_graph(graph, rules, domains, classifiers=None):
    """Traverses a rule graph and classifies a domain collection.

    Terminals are lists of rules without children. So, on a terminal,
    test each rule, breaking on (and saving) the first satisfied.
    Otherwise, test the current node and traverse its children.

    Args:
        graph (list, dict): Rule graph to traverse.
        rules (dict): Rule objects to evaluate on domains.
        domains (list): Domain objects to classify.
        classifiers (list): Current classifiers for a Domain collection.
    """
    if not classifiers:
        classifiers = []

    if isinstance(graph, list):
        for node in graph:
            if isinstance(node, dict):
                classifiers = _traverse_graph(node, rules, domains, classifiers)
                if classifiers:
                    break
            else:
                if rules[node].satisfied_by(domains):
                    rules[node].rename_domains(domains)
                    classifiers.append(node)
                    break
    else:
        for node, children in graph.items():
            if rules[node].satisfied_by(domains):
                rules[node].rename_domains(domains)
                classifiers.append(node)
                classifiers = _traverse_graph(children, rules, domains, classifiers)
                break
    return classifiers


class RuleGraph:
    """Hierarchy of classification rules.

    The RuleGraph is used to classify synthases based on their domains.
    It stores Rule objects, as well as a directed graph controlling the
    order and hierarchy of classification.

    An example synthaser rule graph looks like:
    [
        "Hybrid",
        {"PKS": ["HR-PKS", "PR-PKS", "NR-PKS"]},
        "NRPS"
    ]

    In this example, the "Hybrid" rule is evaluated first. If unsuccessful,
    the "PKS" rule is evaluated. If this is successful, synthaser recurses
    into child rules, in which case the "HR-PKS", "PR-PKS" and "NR-PKS" rules
    can be evaluated, and so on. Each rule name must have a corresponding
    entry in the rules attribute.

    Note that terminal leaves in the graph are placed in lists, whereas
    hierarchies are written as dictionaries of lists. This preserves rule
    order in Python, as well as preventing empty, unnecessary dictionaries
    at every level.

    Attributes:
        rules (dict): Collection of synthaser rules.
        graph (dict): Hierarchy of synthaser rules for classification.
    """

    def __init__(self, rules=None, graph=None):
        self.rules = rules if rules else {}
        self.graph = graph if graph else []

    @classmethod
    def from_dict(cls, d):
        return cls(
            rules={rule["name"]: Rule(**rule) for rule in d["rules"]},
            graph=d["graph"]
        )

    def to_dict(self):
        return {
            "rules": [rule.to_dict() for rule in self.rules],
            "graph": self.graph
        }

    def classify(self, domains):
        return _traverse_graph(self.graph, self.rules, domains)


class Rule:
    """A classification rule.

    Attributes:
        name (str): Name given to proteins satisfying this rule.
        domains (list): Domain types required to satisfy rule.
        filters (dict): Specific CDD families for each domain type.
        evaluator (str): Evaluatable rule satisfaction statement.
    """

    def __init__(self, name=None, rename=None, domains=None, filters=None, evaluator=None):
        self.name = name if name else ""
        self.rename = rename if rename else {}
        self.domains = domains if domains else []
        self.filters = filters if filters else {}
        self.evaluator = evaluator if evaluator else ""

    def to_dict(self):
        return {
            "name": self.name,
            "rename": self.rename,
            "domains": self.domains,
            "filters": self.filters,
            "evaluator": self.evaluator
        }

    def evaluate(self, conditions):
        """Evaluates the rules evaluator string given evaluated conditions.

        Iterates backwards to avoid bad substitutions in larger (>=10) indices.
        e.g. "0 and 1 and ... and 13" --> "False and True and ... and True3"
        """
        evaluator = self.evaluator
        for idx, condition in reversed(list(enumerate(conditions))):
            evaluator = evaluator.replace(str(idx), str(condition))
        return eval(evaluator)

    def rename_domains(self, domains):
        """Renames domain types if substitutions are specified in the rule.

        If the rename dictionary is empty, no action is taken.
        """
        if not self.rename:
            return
        for domain in domains:
            if domain.type in self.rename:
                domain.type = self.rename[domain]

    def valid_family(self, domain):
        """Checks a given domain matches a specified CDD family in the rule.

        If no families have been specified for the given domain type, this
        function will return True (i.e. any family of the type is accepted).
        """
        if domain.type in self.filters:
            return domain.domain in self.filters[domain.type]
        return True

    def satisfied_by(self, domains):
        """Evaluates this rule against a collection of domains.

        Checks that:
        1) required domain types are represented in the supplied domains, and
        2) domains are of the desired CDD families, if any are specified.

        Placeholders in the evaluator string are then replaced by their
        respective booleans, and evaluated.

        Once a domain in the supplied domains has matched one in the rule, it
        cannot be matched to another in the rule. This enables rules based on
        counts of domains (e.g. multi-modular PKS w/ 2 KS domains).
        """
        LOG.debug("Evaluating %s against %s", self.name, [d.type for d in domains])
        seen = []
        conditions = []
        for rule_domain in self.domains:
            match = False
            for domain in domains:
                if domain in seen:
                    continue
                if domain.type == rule_domain and self.valid_family(domain):
                    seen.append(domain)
                    match = True
                    break
            conditions.append(match)
        return self.evaluate(conditions)


def classify(synthases, rule_path):
    with open(rule_path) as fp:
        d = json.load(fp)
        rg = RuleGraph.from_dict(d)
    etc = []
    for synthase in synthases:
        classification = rg.classify(synthase.domains)
        setattr(synthase, "classification", classification)
        etc.append(classification)
    return etc



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

    module = "NRPS" if synthase.type == "NRPS" else "PKS"
    replace = {"ACP": "T", "TR": "R"}

    for domain in synthase.domains:
        # Determine if we're in the PKS or NRPS module
        if domain.type in ("A", "C"):
            module = "NRPS"
        elif domain.type in ("KS", "AT"):
            module = "PKS"

        # Only change domain names if we're in the NRPS module
        if module == "NRPS" and domain.type in replace:
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


def classify_old(synthase):
    """Classify a Synthase.

    First, assign a broad biosynthetic type to the Synthase (PKS, NRPS, Hybrid).
    Then, further classification is performed on Type I PKSs.
    """
    try:
        synthase.type = assign_broad_type(synthase.domains)
    except ValueError:
        return

    if synthase.type == "Type I PKS":
        synthase.subtype = assign_T1PKS_subtype(synthase.domain_types)
    else:
        if synthase.type in ("Hybrid", "NRPS"):
            rename_NRPS_domains(synthase)
        synthase.subtype = synthase.type

"""Create/manipulate domain rules."""

import json
import gzip
import shutil
import tempfile
import logging

from ftplib import FTP
from pathlib import Path
from collections import defaultdict

from synthaser.rpsblast import download as ftp_download


LOG = logging.getLogger(__name__)


def group_domains(domains):
    """Group domain families by their types."""
    d = defaultdict(list)
    for family in domains.values():
        d[family["type"]].append(family["name"])
    return d


def summary(domain_file, rules, families):
    """Generates summaries of the specified domain rules file.

    If no rules or families are specified, the entire file will be summarised.
    Otherwise, the full JSON entries of matching domain rules will be printed.
    """
    # Load in domain rules
    with open(domain_file) as fp:
        df = json.load(fp)

    # If nothing specified, print out summary of entire file
    if not (rules or families):
        types = group_domains(df)
        offset = max(len(key) for key in types)
        for key, value in types.items():
            print(f"{key.rjust(offset)}: {', '.join(value)}")
        return

    # Otherwise, print out records matching rules/families
    for value in df.values():
        if (
            value["type"] in rules
            or value["accession"] in families
            or value["name"] in families
        ):
            print(json.dumps(value, indent=2))


def print_domains(d, domains):
    for item in d["classes"]:
        if item["type"] in domains:
            x = json.dumps(item, indent=2)
            print(x)


def delete_domains(d, domains):
    d["classes"] = [
        domain
        for domain in d["classes"]
        if domain["type"] not in domains
    ]


def download_cdd_families(folder):
    """Downloads CDD family information files from NCBI FTP to a folder."""
    folder = Path(folder)

    if not folder.exists() or not folder.is_dir():
        folder.mkdir()

    cddid_gz = folder / "cddid_all.tbl.gz"
    cddid = cddid_gz.with_suffix("")

    bitscore = folder / "bitscore_specific.txt"
    families = folder / "family_superfamily_links"

    LOG.info("Connecting to NCBI FTP")
    with FTP("ftp.ncbi.nih.gov") as ftp:
        ftp.login()
        ftp.cwd("pub/mmdb/cdd")

        for path in [cddid_gz, bitscore, families]:
            ftp_download(ftp, path.name, path)

        LOG.info("Extracting cddid.tbl.gz")
        with gzip.open(cddid_gz, "rb") as f_in, cddid.open("wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        LOG.info("Deleting cddid.tbl.gz")
        cddid_gz.unlink()


def parse_cdd_families(folder):
    """Generates CDD domain family information database.

    Parses three files downloaded from the NCBI FTP:
    1. cddid_all.tbl: Family information, PSSM length
    2. bitscore_specific.txt: Bitscore thresholds for specific hits
    3. family_superfamily_links: Family-Superfamily relationships

    Then dumps to JSON. File can be used when updating the domain rule
    JSON file with this module.
    """
    folder = Path(folder)

    cddid = folder / "cddid_all.tbl"
    bitscore = folder / "bitscore_specific.txt"
    families = folder / "family_superfamily_links"

    d = {}
    LOG.info("Building CDD family file...")
    with cddid.open() as fp:
        parse_cddid(fp, d)
    with bitscore.open() as fp:
        parse_bitscores(fp, d)
    with families.open() as fp:
        parse_families(fp, d)
    return d


def parse_cddid(fp, d):
    """Parses cddid_all.tbl."""
    for line in fp:
        if not line or line.isspace():
            continue
        pssm, accession, name, _, length = line.strip().split("\t")
        d[accession] = dict(
            pssm=int(pssm),
            accession=accession,
            name=name,
            length=int(length),
        )


def parse_bitscores(fp, d):
    """Parses bitscore_specific.txt."""
    for line in fp:
        if not line or line.isspace():
            continue
        _, accession, bitscore = line.strip().split("\t")
        d[accession]["bitscore"] = float(bitscore)


def parse_families(fp, d):
    """Parses family_superfamily_links."""
    for line in fp:
        if not line or line.isspace():
            continue
        family, _, superfamily, _ = line.strip().split("\t")
        d[family]["superfamily"] = superfamily


def parse_skeleton(fp):
    """Parse domain rule skeleton file."""
    rules = {}
    for entry in fp.read().split(">"):
        if not entry:
            continue
        header, *families = entry.strip().split("\n")
        for family in families:
            rules[family] = dict(type=header)
    return rules


def update_rule(old, new):
    """Merge families in two rules"""
    old["name"] = new["name"]

    # If rule has no families, set to new
    if "families" not in old:
        old["families"] = new["families"]
        return

    # Otherwise, check for matching rules and update each
    old["families"].update(new["families"])


def update_domains(
    rules,
    domains,
    type=None,
    families=None,
    file=None,
):
    """Adds new, or updates existing domain rules."""
    if file:
        LOG.info("Updating rules from skeleton file: %s", file.name)
        for key, value in parse_skeleton(file).items():
            info = lookup(domains, key)
            value.update(info)
            accession = value["accession"]
            rules[accession] = value
    else:
        if not (type and families):
            raise ValueError("Expected type and families")

        LOG.info("Updating %s rule", type)
        for family in families:
            entry = lookup(domains, family)

            if not entry:
                LOG.info("Failed to find %s", family)
                continue

            LOG.info("  Adding %s", family)
            entry["type"] = type
            accession = entry["accession"]
            rules[accession] = entry


def lookup(d, key):
    """Finds a family entry in the CDD JSON file."""
    if key in d:
        return d[key]
    for entry in d.values():
        if key in (entry["name"], entry["pssm"]):
            return entry


def find_accession(d, key):
    """Lookup CDD family accession given a name or accession."""
    return lookup(d, key)["accession"]


def convert_domains(d, names=None, accessions=None):
    """Converts a list of CDD family names to accessions, or vice versa."""
    if not (names or accessions):
        raise ValueError("Expected names or accessions")
    return [
        find_accession(d, item)
        for item in (names or accessions)
    ]


def convert(domain_file, names=None, accessions=None):
    """Converts a list of CDD family names to accessions, or vice versa."""
    d = json.load(domain_file)
    results = convert_domains(d, names=names, accessions=accessions)
    for a, b in zip(names or accessions, results):
        if b:
            print(b)
        else:
            print(f"Failed to find: {a}")


def update(
    rule_file,
    domain_file,
    type=None,
    families=None,
    file=None
):
    """Updates rules in a domain rules file.

    Can either provide a single rule by specifying the type, name and CDD families,
    or the path to a skeleton file for multiple rules.

    If a rule is new (does not exist in existing file), it will be added.
    If a rule exists, but there are new families or a new name, it will be updated.
    """
    LOG.info("Starting domain rule updater module")

    rule_path = Path(rule_file)
    if rule_path.exists():
        LOG.info("Loading domain rules: %s", rule_path.name)
        with rule_path.open() as fp:
            df = json.load(fp)
    else:
        LOG.info("Creating new domain rules file")
        df = {}

    LOG.info("Loading CDD families: %s", domain_file.name)
    db = json.load(domain_file)

    update_domains(
        df,
        db,
        type=type,
        families=families,
        file=file,
    )

    LOG.info("Writing updated rules: %s", rule_path.name)
    with rule_path.open("w") as fp:
        json.dump(df, fp, indent=2)

    LOG.info("Done!")


def remove(rule_file, rules, families):
    """Removes rules or CDD families from a domain rules file.

    If multiple rules are specified, their entire entries will be removed from the file.
    If a single rule, but multiple families are specified, those families will be
    removed from the specific rule.
    """
    LOG.info("Starting domain removal module")
    rule_path = Path(rule_file)

    LOG.info("Loading domain rules: %s", rule_path.name)
    with rule_path.open() as fp:
        df = json.load(fp)

    if rules:
        LOG.info("Deleting rule types: %s", rules)
    if families:
        LOG.info("Deleting families: %s", families)

    for key in list(df):
        entry = df[key]
        if (
            entry["type"] in rules
            or entry["accession"] in families
            or entry["name"] in families
        ):
            LOG.info(" %s, %s [%s]", entry["type"], entry["name"], entry["accession"])
            df.pop(key)

    LOG.info("Writing updated rules: %s", rule_path.name)
    with rule_path.open("w") as fp:
        json.dump(df, fp, indent=2)

    LOG.info("Done!")


def download(domain_file, folder=None, indent=2):
    """Download CDD family information."""
    LOG.info("Starting CDD family download module")
    if folder:
        download_cdd_families(folder)
        d = parse_cdd_families(folder)
    else:
        LOG.info("Using a temporary directory.")
        with tempfile.TemporaryDirectory() as tmpdir:
            download_cdd_families(tmpdir)
            d = parse_cdd_families(tmpdir)
    LOG.info("Writing: %s" % domain_file)
    with open(domain_file, "w") as fp:
        json.dump(d, fp, indent=indent)
    LOG.info("Done!")

"""Extract CDS sequences from GenBank files for synthaser analysis."""

import logging

from pathlib import Path

from Bio import SeqIO


LOG = logging.getLogger(__name__)


def get_feature_type(feature):
    """Tests if an antiSMASH feature is a NRPS or PKS.

    Looks for a /NRPS_PKS="type:..." field containing either 'PKS' or 'NRPS'.
    This function returns True if a matching qualifier is found.
    """
    try:
        for value in feature.qualifiers["NRPS_PKS"]:
            if not value.startswith("type:"):
                continue
            if "PKS" in value:
                return "PKS"
            if "NRPS" in value:
                return "NRPS"
    except KeyError:
        return None


def get_NRPS_PKS(features):
    pks, nrps = [], []
    for feature in features:
        ftype = get_feature_type(feature)
        if ftype == "PKS":
            pks.append(feature)
        elif ftype == "NRPS":
            nrps.append(ftype)
    return pks, nrps


def convert(path, antismash=False):
    """Convert a GenBank file to FASTA.

    Arguments:
        path (str): GenBank file path
        output (str): Output file path
        antismash (bool): Only save PKS/NRPS sequences
    Returns:
        records (list): genome2json Feature objects of parsed records
    """
    LOG.info("Parsing GenBank file: %s", path)

    path = Path(path)

    with path.open() as fp:
        features = [record for record in SeqIO.parse(fp, "genbank-cds")]

    if antismash:
        LOG.info("Finding antiSMASH PKS and NRPS features")
        pks, nrps = get_PKS_NRPS(features)
        for fts, text in [(pks, 'pks'), (nrps, 'nrps')]:
            with path.with_name(f"{path.name}_{text}.fa").open("w") as fp:
                LOG.info("Writing %s: %s", text.upper(), fp.name)
                SeqIO.write(fts, fp, "fasta")
    else:
        with path.with_suffix(".fa").open("w") as fp:
            LOG.info("Writing FASTA: %s", fp.name)
            SeqIO.write(features, fp, "fasta")

    return features

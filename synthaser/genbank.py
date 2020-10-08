"""Extract CDS sequences from GenBank files for synthaser analysis."""

import logging

from g2j import genbank


LOG = logging.getLogger(__name__)


def fasta(features):
    """Builds a FASTA string from g2j Feature objects."""
    return "\n".join(
        ">{}\n{}".format(
            f.qualifiers.get("locus_tag")
            or f.qualifiers.get("protein_id"),
            f.qualifiers.get("translation")
        )
        for f in features
    )


def feature_is_NRPS_PKS(feature):
    """Tests if an antiSMASH feature is a NRPS or PKS.

    Looks for a /NRPS_PKS="type:..." field containing either 'PKS' or 'NRPS'.
    This function returns True if a matching qualifier is found.
    """
    if "NRPS_PKS" not in feature.qualifiers:
        return False
    return any(
        "PKS" in value or "NRPS" in value
        for value in feature.qualifiers["NRPS_PKS"]
        if value.startswith("type:")
    )


def parse(handle, antismash=False):
    """Parses a GenBank file handle for sequence features.

    If antismash is True, any features not containing a /NRPS_PKS="type:"
    field containing 'PKS' or 'NRPS' is discarded.
    """
    organism = genbank.parse(
        handle,
        feature_types=["CDS"],
        save_scaffold_sequence=False
    )

    features = [
        feature
        for scaffold in organism.scaffolds
        for feature in scaffold.features
    ]

    if antismash:
        # Filter out any non-PKS/NRPS in an antiSMASH file
        LOG.info("antiSMASH file specified; finding PKS/NRPS")
        features = [
            feature
            for feature in features
            if feature_is_NRPS_PKS(feature)
        ]

    return features


def convert(path, output=None, antismash=False):
    """Convert a GenBank file to FASTA.

    Arguments:
        path (str): GenBank file path
        output (str): Output file path
        antismash (bool): Only save PKS/NRPS sequences
    Returns:
        records (list): genome2json Feature objects of parsed records
    """
    LOG.info("Parsing GenBank file: %s", path)
    with open(path) as fp:
        records = parse(fp, antismash=antismash)

    result = fasta(records)

    if output:
        LOG.info("Writing FASTA file: %s", output)
        with open(output, "w") as fp:
            fp.write(result)

    return result

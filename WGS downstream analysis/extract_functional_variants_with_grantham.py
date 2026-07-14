#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import argparse
import pandas as pd

# ----------------------------------------------------------
# Complete Grantham matrix (20 amino acids, one-letter code)
# ----------------------------------------------------------
GRANTHAM = {
    ('A','R'):112, ('A','N'):111, ('A','D'):126, ('A','C'):195, ('A','Q'):91,  ('A','E'):107, ('A','G'):60,  ('A','H'):86,  ('A','I'):94,  ('A','L'):96,
    ('A','K'):106, ('A','M'):84,  ('A','F'):113, ('A','P'):27,  ('A','S'):99,  ('A','T'):58,  ('A','W'):148, ('A','Y'):112, ('A','V'):64,

    ('R','N'):86,  ('R','D'):96,  ('R','C'):180, ('R','Q'):43,  ('R','E'):54,  ('R','G'):125, ('R','H'):29,  ('R','I'):97,  ('R','L'):102, ('R','K'):26,
    ('R','M'):91,  ('R','F'):97,  ('R','P'):103, ('R','S'):110, ('R','T'):71,  ('R','W'):101, ('R','Y'):77,  ('R','V'):96,

    ('N','D'):23,  ('N','C'):139, ('N','Q'):46,  ('N','E'):42,  ('N','G'):80,  ('N','H'):68,  ('N','I'):149, ('N','L'):153, ('N','K'):94,  ('N','M'):142,
    ('N','F'):158, ('N','P'):91,  ('N','S'):46,  ('N','T'):65,  ('N','W'):174, ('N','Y'):143, ('N','V'):133,

    ('D','C'):154, ('D','Q'):61,  ('D','E'):45,  ('D','G'):94,  ('D','H'):81,  ('D','I'):168, ('D','L'):172, ('D','K'):101, ('D','M'):160, ('D','F'):177,
    ('D','P'):108, ('D','S'):65,  ('D','T'):85,  ('D','W'):181, ('D','Y'):160, ('D','V'):152,

    ('C','Q'):154, ('C','E'):170, ('C','G'):159, ('C','H'):174, ('C','I'):198, ('C','L'):198, ('C','K'):202, ('C','M'):196, ('C','F'):205, ('C','P'):169,
    ('C','S'):112, ('C','T'):149, ('C','W'):215, ('C','Y'):194, ('C','V'):192,

    ('Q','E'):29,  ('Q','G'):87,  ('Q','H'):24,  ('Q','I'):109, ('Q','L'):113, ('Q','K'):53,  ('Q','M'):101, ('Q','F'):116, ('Q','P'):76,  ('Q','S'):68,
    ('Q','T'):42,  ('Q','W'):130, ('Q','Y'):99,  ('Q','V'):96,

    ('E','G'):98,  ('E','H'):40,  ('E','I'):134, ('E','L'):138, ('E','K'):56,  ('E','M'):126, ('E','F'):140, ('E','P'):93,  ('E','S'):80,  ('E','T'):65,
    ('E','W'):152, ('E','Y'):122, ('E','V'):121,

    ('G','H'):98,  ('G','I'):135, ('G','L'):138, ('G','K'):127, ('G','M'):127, ('G','F'):153, ('G','P'):42,  ('G','S'):56,  ('G','T'):59,  ('G','W'):184,
    ('G','Y'):147, ('G','V'):109,

    ('H','I'):94,  ('H','L'):99,  ('H','K'):32,  ('H','M'):87,  ('H','F'):100, ('H','P'):77,  ('H','S'):89,  ('H','T'):47,  ('H','W'):115, ('H','Y'):83,
    ('H','V'):84,

    ('I','L'):5,   ('I','K'):102, ('I','M'):10,  ('I','F'):21,  ('I','P'):95,  ('I','S'):142, ('I','T'):89,  ('I','W'):61,  ('I','Y'):33,  ('I','V'):29,

    ('L','K'):107, ('L','M'):15,  ('L','F'):22,  ('L','P'):98,  ('L','S'):145, ('L','T'):92,  ('L','W'):61,  ('L','Y'):36,  ('L','V'):32,

    ('K','M'):95,  ('K','F'):102, ('K','P'):103, ('K','S'):121, ('K','T'):78,  ('K','W'):110, ('K','Y'):85,  ('K','V'):97,

    ('M','F'):28,  ('M','P'):87,  ('M','S'):135, ('M','T'):81,  ('M','W'):67,  ('M','Y'):36,  ('M','V'):21,

    ('F','P'):114, ('F','S'):155, ('F','T'):103, ('F','W'):40,  ('F','Y'):22,  ('F','V'):50,

    ('P','S'):74,  ('P','T'):38,  ('P','W'):147, ('P','Y'):110, ('P','V'):68,

    ('S','T'):58,  ('S','W'):177, ('S','Y'):144, ('S','V'):124,

    ('T','W'):128, ('T','Y'):92,  ('T','V'):69,

    ('W','Y'):37,  ('W','V'):88,

    ('Y','V'):55
}

AA3_TO_AA1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
    'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
    'Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
    'Ter':'*','Stop':'*','Sec':'U','Pyl':'O'
}

LOF_EFFECTS = {
    'stop_gained',
    'frameshift_variant',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'start_lost',
    'stop_lost'
}

SYN_EFFECTS = {
    'synonymous_variant',
    'start_retained_variant',
    'stop_retained_variant'
}

MISSENSE_EFFECTS = {
    'missense_variant'
}

def get_grantham(a1, a2):
    if a1 is None or a2 is None:
        return None
    if a1 == a2:
        return 0
    if a1 == '*' or a2 == '*':
        return None
    if (a1, a2) in GRANTHAM:
        return GRANTHAM[(a1, a2)]
    if (a2, a1) in GRANTHAM:
        return GRANTHAM[(a2, a1)]
    return None

def parse_hgvs_p(hgvs_p):
    """
    Parse HGVS protein change examples:
      p.Ala123Asp
      p.Gly45Trp
      p.Trp24Ter
      p.Ser77=
    Returns: (aa_ref, aa_alt, aa_pos)
    """
    if hgvs_p is None or hgvs_p == '' or hgvs_p == '.':
        return None, None, None

    # standard amino acid substitution
    m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|Stop|=)', hgvs_p)
    if m:
        ref3, pos, alt3 = m.group(1), m.group(2), m.group(3)
        aa_ref = AA3_TO_AA1.get(ref3, None)
        aa_alt = '=' if alt3 == '=' else AA3_TO_AA1.get(alt3, None)
        return aa_ref, aa_alt, pos

    return None, None, None

def parse_ann_record(ann):
    """
    SnpEff ANN format (standard):
    Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID |
    Transcript_BioType | Rank/Total | HGVS.c | HGVS.p | cDNA.pos / cDNA.length |
    CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
    """
    fields = ann.split('|')
    # pad short fields
    if len(fields) < 16:
        fields += [''] * (16 - len(fields))

    return {
        'allele': fields[0],
        'effect': fields[1],
        'impact': fields[2],
        'gene': fields[3],
        'gene_id': fields[4],
        'feature_type': fields[5],
        'transcript': fields[6],
        'biotype': fields[7],
        'rank_total': fields[8],
        'HGVS.c': fields[9],
        'HGVS.p': fields[10],
        'cDNA_pos_len': fields[11],
        'CDS_pos_len': fields[12],
        'AA_pos_len': fields[13],
        'distance': fields[14],
        'info': fields[15]
    }

def classify_variant(effect, hgvs_p, severe_threshold=150):
    aa_ref, aa_alt, aa_pos = parse_hgvs_p(hgvs_p)
    grantham_score = None
    var_class = None

    if effect in LOF_EFFECTS:
        var_class = 'LoF'

    elif effect in MISSENSE_EFFECTS:
        var_class = 'missense_all'
        if aa_ref is not None and aa_alt is not None and aa_alt != '=':
            grantham_score = get_grantham(aa_ref, aa_alt)
            if grantham_score is not None and grantham_score > severe_threshold:
                var_class = 'missense_severe'

    elif effect in SYN_EFFECTS:
        var_class = 'synonymous'
        if aa_ref is not None and aa_alt == '=':
            grantham_score = 0

    return var_class, aa_ref, aa_alt, aa_pos, grantham_score

def main():
    parser = argparse.ArgumentParser(description="Extract LoF / missense / severe missense (Grantham) / synonymous from SnpEff ANN table.")
    parser.add_argument("-i", "--input", required=True, help="Input ANN TSV generated by bcftools query")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    parser.add_argument("--severe-threshold", type=int, default=150, help="Grantham score threshold for severe missense (default: 150)")
    args = parser.parse_args()

    records = []

    with open(args.input) as f:
        for line_num, line in enumerate(f, start=1):
            line = line.rstrip('\n')
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 5:
                sys.stderr.write(f"[WARN] line {line_num} has <5 columns, skipped\n")
                continue

            chrom, pos, ref, alt, ann_field = parts[0], parts[1], parts[2], parts[3], parts[4]

            if ann_field == '.' or ann_field == '':
                continue

            ann_records = ann_field.split(',')

            for ann in ann_records:
                parsed = parse_ann_record(ann)
                effect = parsed['effect']
                hgvs_p = parsed['HGVS.p']

                var_class, aa_ref, aa_alt, aa_pos, grantham_score = classify_variant(
                    effect,
                    hgvs_p,
                    severe_threshold=args.severe_threshold
                )

                if var_class is None:
                    continue

                records.append({
                    'CHROM': chrom,
                    'POS': pos,
                    'REF': ref,
                    'ALT': alt,
                    'effect': effect,
                    'impact': parsed['impact'],
                    'gene': parsed['gene'],
                    'gene_id': parsed['gene_id'],
                    'feature_type': parsed['feature_type'],
                    'transcript': parsed['transcript'],
                    'biotype': parsed['biotype'],
                    'rank_total': parsed['rank_total'],
                    'HGVS.c': parsed['HGVS.c'],
                    'HGVS.p': parsed['HGVS.p'],
                    'AA_ref': aa_ref,
                    'AA_alt': aa_alt,
                    'AA_pos': aa_pos,
                    'Grantham': grantham_score,
                    'class': var_class
                })

    df = pd.DataFrame(records)

    if df.empty:
        sys.stderr.write("[WARN] No functional variants extracted.\n")
        df.to_csv(args.output, sep='\t', index=False)
        return

    # Optional: remove exact duplicate rows
    df = df.drop_duplicates()

    # Sort
    try:
        df['POS'] = df['POS'].astype(int)
    except:
        pass

    df = df.sort_values(by=['CHROM', 'POS', 'gene', 'transcript', 'effect'])

    df.to_csv(args.output, sep='\t', index=False)

    # Summary to stderr
    sys.stderr.write("=== Summary ===\n")
    sys.stderr.write(df['class'].value_counts(dropna=False).to_string() + "\n")
    severe_n = (df['class'] == 'missense_severe').sum()
    sys.stderr.write(f"Severe missense (Grantham > {args.severe_threshold}): {severe_n}\n")
    sys.stderr.write(f"Output written to: {args.output}\n")

if __name__ == "__main__":
    main()

"""
IO.vcf

Serialization utilities for writing DotPlotLineSegment model objects
to VCF format using pysam.
"""

import pysam

def write_vcf(SVs, output_path="output.vcf", sample='SAMPLE'):
    header = pysam.VariantHeader()
    header.add_meta('INFO', items=[('ID', 'SVTYPE'), ('Number', '1'), ('Type', 'String'), ('Description', 'Type of structural variant')])
    header.add_meta('INFO', items=[('ID', 'END'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'End position of the SV')])
    header.add_meta('INFO', items=[('ID', 'SVLEN'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Length of the SV')])
    header.add_meta('INFO', items=[('ID', 'THETA'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Angle between segments in degrees')])
    header.add_meta('INFO', items=[('ID', 'IMPRECISE'), ('Number', '0'), ('Type', 'Flag'), ('Description', 'Imprecise structural variation')])
    header.add_meta('INFO', items=[('ID', 'DUP_TANDEM'), ('Number', '0'), ('Type', 'Flag'), ('Description', 'Tandem duplication')])
    header.add_meta('INFO', items=[('ID', 'CN'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Estimated copy number in query')])
    header.add_meta('INFO', items=[('ID', 'CNL'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Estimated copy number in reference')])
    header.add_meta('INFO', items=[('ID', 'EVENT'),('Number', '1'),('Type', 'String'),('Description', 'ID of the complex event this variant is part of')])
    
    header.add_sample(sample)

    contigs = sorted(set(SV.chrom for SV in SVs)) # get each unique contig
    for contig in contigs:
        header.contigs.add(contig) ## add them to the vcf header

    vcf_out = pysam.VariantFile(output_path, 'w', header=header)

    for i, SV in enumerate(SVs): #########################
        if SV.theta == 45 and not (SV.range_partners or SV.domain_partners):
            continue
        chrom = SV.chrom

        if SV.theta < -90 or SV.theta > 90:
            pos = SV.reference_end
            end = SV.reference_start
        else:
            pos = SV.reference_start
            end = SV.reference_end
    
        length = end - pos if SV.sv_type != "INS" else SV.query_end - SV.query_start
        ref_base = 'N'
        alt = f"<{SV.sv_type if SV.sv_type != 'DUP_REF' else 'DEL'}>"

        info_dict = {
            'SVTYPE': SV.sv_type if SV.sv_type != 'DUP_REF' else 'DEL',
            'END': end,
            'SVLEN': length,
            'EVENT': SV.event_ID,
            'THETA': SV.theta,
            'IMPRECISE': True
        }

        # Annotate tandem duplication flags
        if SV.sv_type in ('DUP', 'DUP_REF'):
            info_dict['DUP_TANDEM'] = True
            info_dict['CN'] = 2 if SV.sv_type == 'DUP' else 1
            info_dict['CNL'] = 1 if SV.sv_type == 'DUP' else 2

        record = vcf_out.new_record(
            contig=chrom,
            start=pos - 1,
            stop=end,
            alleles=(ref_base, alt),
            id=f'sv_{i}',
            info=info_dict
        )
        vcf_out.write(record)

    vcf_out.close()
#! /usr/bin/env python3
import sys
import argparse
import vcf

def parse_ped(ped_file):
    samples = {}  # 'sample_id': {'family': 'fam_id', 'parents': ['sample_id', 'sample_id']}
    for line in ped_file:
        ped_data = line.strip().split()
        family, sample, father, mother, sex, phenotype = ped_data

        # Create samples
        if sample not in samples:
            samples[sample] = {'family': family, 'parents': [], 'children': []}
        if father != '0' and father not in samples:
            samples[father] = {'family': family, 'parents': [], 'children': []}
        if mother != '0' and mother not in samples:
            samples[mother] = {'family': family, 'parents': [], 'children': []}

        # Save sample relations
        if father != '0':
            samples[sample]['parents'].append(father)
            samples[father]['children'].append(sample)
        if mother != '0':
            samples[sample]['parents'].append(mother)
            samples[mother]['children'].append(sample)

    families = {}
    for sample in samples:
        if len(samples[sample]['parents']) == 2:
            families[sample] = samples[sample]['parents']

    return samples, families


def parse_vcf(vcf_file):  # returns list with genotypes
    with open(vcf_file, 'r') as vcf_input_file:
        vcf_reader = vcf.Reader(vcf_input_file)

        if 'fileformat' not in vcf_reader.metadata: # Check if true VCF file
            sys.exit("Input file {} is not a correct VCF file. "\
                "Field \"fileformat\" was not detected in VCF".format(args.inputfile)
                )
        if len(vcf_reader.samples) > 1: # Check is VCF is not a multisample VCF
            sys.exit("Single sample VCF support only. Input file {} is a multisample VCF".format(args.inputfile)) 

        sampleid = vcf_reader.samples[0]
        snv_list = []
        for record in vcf_reader:
            chrom = record.CHROM
            pos = record.POS

            if 'DP' in record.genotype(sampleid).data._asdict():
                dp = record.genotype(sampleid)['DP']
            else: # Skip variant position if DP is missing
                continue

            if dp >= args.mindepth and chrom is not "Y" :
                if record.genotype(sampleid).phased: # sort genotype in unphased state.
                    gt = record.genotype(sampleid)['GT'].split("|")
                    gt.sort()
                    gt = "/".join(gt)
                else:
                    gt = record.genotype(sampleid)['GT']
                if len(record.ALT) > 1: # skip multiallelic position
                   continue
                allele = [str(record.REF[0]), str(record.ALT[0])]
                genotype_call = record.genotype(sampleid).is_variant
                if genotype_call is not None: # True = ref/alt or alt/alt genotype, False = ref/ref genotype, None = missing genotype (i.e. ./.). Skip missing genotypes in statement.
                    genotype = [] 
                    for item in gt.split("/"):
                        genotype.append(allele[int(item)])
                    allele_genotype = "/".join(genotype)
                    snv_list.append(["{}_{}".format(chrom, pos), gt])
    return snv_list


def make_upd(families, samples):
    vcfs = {}
    for vcf in args.input_files:
        sampleid = vcf.split(args.suffix)[0].split("/")[-1].strip("_dedup.realigned")
        if sampleid not in vcfs:
            vcfs[sampleid] = vcf  
    for sample in families:
        family = samples[sample]['family']

        if sample not in vcfs:  # PED file includes all samples of entire run, not only project. Skip if child VCF is not present in current analysis.
            continue
    
        if args.sample_id not in sample:
            continue

        child = parse_vcf(vcfs[sample])
        father = dict(parse_vcf(vcfs[families[sample][0]]))  ## father always first item in families dict
        mother = dict(parse_vcf(vcfs[families[sample][1]]))  ## mother always second item in families dict

        output_file = open("{}_{}.igv".format(args.run_id, family), 'w')
        output_file.write("#track type=igv name=Pat/Mat_segments color=204,204,0 altColor=0,100,224 graphType=bar windowingFunction=none maxHeightPixels=50 viewLimits=-1,1\n")
        chromosome = "1" 
        start = 0
        for variant in child:
            if variant[0] in father and variant[0] in mother:  # position is called in both father and mother
                position, child_geno, father_geno, mother_geno = variant[0], variant[1], father[variant[0]], mother[variant[0]]
                label = '' 
                score = ''
                chrom = str(position.split("_")[0])
                if chrom != chromosome:
                    start = 0  
                    chromosome = chrom

                snv_pos = int(position.split("_")[1])
                start = start
                end = snv_pos

                genotype_conversion = {"0/0":"homref","0/1":"het", "1/1":"homvar"}

                genotype_score = {
                    "homref_homvar_homref": ["patIso", 1],
                    "homref_homvar_homvar": ["matIso", -1],
                    "homref_homvar_het": ["normal", -0.2],
                    "homvar_homref_homvar": ["patIso", 1],
                    "homvar_homref_homref": ["matIso", -1],
                    "homvar_homref_het": ["normal", 0.2],
                    "homvar_het_homref": ["matIso", -1],
                    "het_homvar_homref": ["patIso", 1],
                    "het_homref_het": ["patHet", 1],
                    "homref_het_het": ["matHet", -1],
                    "het_homref_homvar": ["patIso", 1], 
                    "homref_het_homvat": ["matIso", -1]
                    }

                if father_geno in genotype_conversion and mother_geno in genotype_conversion and child_geno in genotype_conversion:
                    genotype = "{}_{}_{}".format(genotype_conversion[father_geno], genotype_conversion[mother_geno],genotype_conversion[child_geno])
                    if genotype in genotype_score:
                        label = genotype_score[genotype][0]
                        score = genotype_score[genotype][1]

                binsize = 0
            
                if int(end - start) > args.maxlocus:
                    start = end - args.maxlocus

                if label: 
                    tag = "Name={label};Position={chromosome}:{postion};Child={child_geno};Father={father_geno},Mother={mother_geno}".format(
                        label=label,
                        chromosome=chrom,
                        postion=snv_pos, 
                        child_geno=child_geno,
                        father_geno=father_geno,
                        mother_geno=mother_geno
                        )

                    output_file.write("{chrom}\t{start}\t{end}\t{tag}\t{score}\n".format(
                        chrom=chrom,
                        start=start-binsize,
                        end=end+binsize,
                        tag=tag,
                        score=score
                        ))

                    start = snv_pos
                   

        output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('ped_file', type=argparse.FileType('r'), help='PED file')
    parser.add_argument('run_id', help='run ID prefix')
    parser.add_argument('sample_id', help='sample id to be processed')
    parser.add_argument('input_files', nargs='+', help='input files (space separated)')
    parser.add_argument('--mindepth', default= 15, type=int, help='Threshold for minimum depth (DP) of SNV (default = 15)')
    parser.add_argument('--suffix', default= ".vcf", type=str, help='suffix of VCF file to be searched (default = .vcf)')
    parser.add_argument('--maxlocus', default= 1000000, type=int, help='maximum size of locus to be printed. This reduces large blocks in regions with low informativity (default = 1000000)')
    args = parser.parse_args()

    samples, families = parse_ped(args.ped_file)
    if args.sample_id in families:
        if not [sampleids for sampleids in args.input_files if args.sample_id in sampleids]: ## If sample_id has VCF in --input_files
            sys.exit("ERROR: VCF of sample {} is missing or misspelled in sample_id argument.".format(args.sample_id))

        parent_missing = []
        for parent in families[args.sample_id]: # check if both parents (based on pedigree file) are present in --input_files
            if not [sampleids for sampleids in args.input_files if parent in sampleids]:
                parent_missing.append(parent)
        if parent_missing:
            sys.exit("ERROR: VCF of parent(s) {} is missing in --input_files, or the pedigree file is incorrect.".format(parent_missing))

    make_upd(families, samples)

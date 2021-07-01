#! /usr/bin/env python3
import subprocess
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
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        sampleid = vcf_reader.samples[0]
        snv_list = []
        for record in vcf_reader:
            chrom = record.CHROM
            pos = record.POS
            dp = record.genotype(sampleid)['DP']
            if dp >= args.mindepth and chrom is not "Y" :
                if record.genotype(sampleid).phased: # sort genotype in unphased state.
                    gt = record.genotype(sampleid)['GT'].split("|")
                    gt.sort()
                    gt = "/".join(gt)
                else:
                    gt = record.genotype(sampleid)['GT']
                if len(record.ALT) > 1: # skip multiallelic position
                   continue
                allele = [str(record.REF[0]),str(record.ALT[0])]
                variant_call = record.genotype(sampleid).is_variant
                if variant_call is not None:
                    genotype = [] 
                    for item in gt.split("/"):
                        genotype.append(allele[int(item)])
                    allele_genotype = "/".join(genotype)
                    snv_list.append(["{}_{}".format(chrom, pos), gt])
    return snv_list


def make_upd(families, samples):
    vcf_files = subprocess.getoutput("find -L {input} -type f -iname \"*{suffix}\"".format(input=args.inputfolder, suffix=args.suffix)).split()  ## Make input parameter in NF?
    print(vcf_files)
    colors = ["204,204,0", "0,100,224","192,192,192"]
    vcfs = {}
    for vcf in vcf_files:
        sampleid = vcf.split(args.suffix)[0].split("/")[-1].strip("_dedup.realigned")
        if sampleid not in vcfs:
            vcfs[sampleid] = vcf  
    for sample in families:  ## This assumes that all VCF of all family members are present! Otherwise hard crash....fix this
        family = samples[sample]['family']
        if sample not in vcfs:
            continue
        child = parse_vcf(vcfs[sample])
        father = dict(parse_vcf(vcfs[families[sample][0]]))  ## father always first item
        mother = dict(parse_vcf(vcfs[families[sample][1]]))  ## mother always second item

        output_file = open("{}_{}.igv".format(args.outputfile, family), 'w')
        output_file.write("track type=igv name=UPD_track color=204,204,0 altColor=0,100,224 graphType=bar windowingFunction=none maxHeightPixels=50 viewLimits=-1,1")
        chromosome = "1" 
        start = 0
        for variant in child:
            if variant[0] in father and variant[0] in mother:  # position is called in both father and mother
                position, child_geno, father_geno, mother_geno = variant[0], variant[1], father[variant[0]], mother[variant[0]]
                label = '' 
                rgb = ''
                chrom = str(position.split("_")[0])
                if chrom != chromosome:
                    start = 0  
                    chromosome = chrom

                snv_pos = int(position.split("_")[1])
                start = start
                end = snv_pos

                ## Make dictionary here? 
                if father_geno == "0/0" and mother_geno == "1/1":
                    if child_geno == "0/0":
                        label = "patIso" 
                        rbg = colors[0]
                        score = 1
                    elif child_geno == "1/1":
                        label = "matIso"
                        rbg = colors[1] 
                        score = -1
                    elif child_geno == "0/1":
                        label = "normal"
                        rbg = colors[2]
                        score = 0
                elif father_geno == "1/1" and mother_geno == "0/0":
                    if child_geno == "1/1":
                        label = "patIso"
                        rbg = colors[0]
                        score = 1
                    elif child_geno == "0/0":
                        label = "matIso"
                        rbg = colors[1]
                        score = -1
                    elif child_geno == "0/1":
                        label = "normal"
                        rbg = colors[2]
                        score = 0
                elif father_geno == "1/1" and mother_geno == "0/1":
                    if child_geno == "0/0":
                        label = "matIso"
                        rbg = colors[1]
                        score = -1
                elif father_geno == "0/1" and mother_geno == "1/1":
                    if child_geno == "0/0":
                        label = "patIso"
                        rbg = colors[0]
                        score = 1
                elif father_geno == "0/1" and mother_geno == "0/0":
                    if child_geno == "0/1":
                       label = "patHet"
                       rbg = colors[0]
                       score = 1
                elif father_geno == "0/0" and mother_geno == "0/1":
                    if child_geno == "0/1":
                       label = "matHet"
                       rbg = colors[1]
                       score = -1
                elif father_geno == "0/1" and mother_geno == "0/0":
                    if child_geno == "1/1":
                       label = "patIso"
                       rbg = colors[0] 
                       score = 1 
                elif father_geno == "0/0" and mother_geno == "0/1":
                    if child_geno == "1/1":
                       label = "matIso"
                       rbg = colors[1]
                       score = -1
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

                    #output_file.write("{chrom}\t{start}\t{end}\t{label}\t{score}\t{strand}\t{thickStart}\t{thickEnd}\t{itemRgb}\t{blockCount}\t{blockSizes}\t{blockStarts}\n".format(
                    #    chrom=chrom,
                    #    start=start-binsize,
                    #    end=end+binsize,
                    #    label=label,
                    #    score=0,
                    #    strand=".",
                    #    thickStart=start-binsize,
                    #    thickEnd=end+binsize,
                    #    itemRgb=rbg,
                    #    blockCount=1,
                    #    blockSizes=end-start+(2*binsize),
                    #    blockStarts=0          
                    #    ))
                    start = snv_pos
                   

        output_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to folder including VCF files')
    parser.add_argument('ped_file', type=argparse.FileType('r'), help='PED file')
    parser.add_argument('outputfile', help='output prefix filename (i.e. runID)')
    #parser.add_argument('--gq_thres', default= 99, type=int, help='Threshold for minimum genotypequality (default = 99)')
    parser.add_argument('--mindepth', default= 15, type=int, help='Threshold for minimum depth (DP) of SNV (default = 15)')
    parser.add_argument('--suffix', default= ".vcf", type=str, help='suffix of VCF file to be searched (default = .vvcf)')
    parser.add_argument('--maxlocus', default= 1000000, type=int, help='maximum size of locus to be printed. This reduces large blocks in regions with low informativity (default = 1000000)')
    args = parser.parse_args()

    samples, families = parse_ped(args.ped_file)
    make_upd(families, samples)

process annotate {
	label "VEP"

	cache "lenient"
        scratch true
	cpus 1

	containerOptions "-B ${params.VEP_dir}/vep_cache:/opt/vep/.vep -B ${params.Gen_pipeline_data}/genomes/${params.genome}:/opt/vep/.genome"

	input:
	tuple path(vcf), path(bed)
	
	output:
	path("${vcf.getSimpleName()}.vep.tsv"), emit: Annotations

	publishDir "${params.Auxiliary_files}/${params.Library_Type}/Annotations", pattern: "${vcf.getSimpleName()}.vep", mode: "copy"
	script:
	def species
	def cache_version
    switch (params.genome) {
        case "hg38":
            species = "homo_sapiens"
	    cache_version=''
            break
        case "mm39":
            species = "mus_musculus"
	    cache_version=''
            break
        case "TAIR10":
            species = "arabidopsis_thaliana"
            cache_version='--cache_version 60'
            break
        case "CHOK1S_HZDv1":
            species = "cricetulus_griseus_chok1gshd"
            break
        case "ASM584v2":
            species = "escherichia_coli_o157_h7_str_sakai_gca_000008865"
	    cache_version='--cache_version 60'
            break
        case "GRCg7b":
            species = "gallus_gallus"
	    cache_version=''
            break
        case "mRatBN7.2":
            species = "rattus_norvegicus"
	    cache_version=''
            break
        default:
            species = "Unknown"
            cache_version=''
    }
	if (params.genome == "hg38")
		"""
	editor="${vcf.simpleName.tokenize('_')[3]}"
	awk '{if (\$4 ~ /^custom/) {print > "custom_regions.bed"} else {print \$4 > "protein.list"}}' ${bed}
	export PERL5LIB=/opt/vep/.vep/Plugins/:\${PERL5LIB:-}
	bcftools sort ${vcf} -Oz -o ${vcf.getSimpleName()}.sorted.vcf.gz
        tabix -p vcf ${vcf.getSimpleName()}.sorted.vcf.gz
	
	if [ -e protein.list ] ; then
		vep -o stdout -i ${vcf.getSimpleName()}.sorted.vcf.gz  --flag_pick_allele_gene --cache --offline -species homo_sapiens --assembly GRCh38 --format vcf --tab --force_overwrite --dir_cache /opt/vep/.vep/ --dir /opt/vep/.vep/Plugins --plugin CADD,/opt/vep/.vep/CADD_GRCh38/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/CADD_GRCh38/InDels.tsv.gz --plugin CONTEXT --variant_class --sift b --polyphen b --nearest gene --gene_phenotype --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype  --pubmed --shift_hgvs 0 --allele_number --buffer_size 10000 --custom ${vcf.getSimpleName()}.sorted.vcf.gz,\$editor,vcf,exact,0,Protospacer,PAM,Nchange | filter_vep --filter "SYMBOL in protein.list" -o "${vcf.getSimpleName()}_Protein.vep"
	else
		touch "${vcf.getSimpleName()}_Protein.vep"
	fi
	if [ -e custom_regions.bed ] ; then
		bcftools view ${vcf.getSimpleName()}.sorted.vcf.gz -R custom_regions.bed | vep -o "${vcf.getSimpleName()}_region.vep" --flag_pick --cache --offline -species homo_sapiens --assembly GRCh38 --format vcf --tab --force_overwrite --dir_cache /opt/vep/.vep/ --dir /opt/vep/.vep/Plugins --plugin CADD,/opt/vep/.vep/CADD_GRCh38/whole_genome_SNVs.tsv.gz,/opt/vep/.vep/CADD_GRCh38/InDels.tsv.gz --plugin CONTEXT --variant_class --sift b --polyphen b --nearest gene --gene_phenotype --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory  --custom ${vcf.getSimpleName()}.sorted.vcf.gz,\$editor,vcf,exact,0,Protospacer,PAM,Nchange
	else
		touch ${vcf.getSimpleName()}_region.vep
	fi
	cat "${vcf.getSimpleName()}_Protein.vep" "${vcf.getSimpleName()}_region.vep"  | awk '/^##/ {next} /^#/ && c++ {next} 1' > tmp.tsv
	awk -v value=\$editor 'BEGIN {OFS="\t"} NR==1 {\$0 = \$0 OFS "editor"} NR>1 {\$0 = \$0 OFS value} 1' tmp.tsv > ${vcf.getSimpleName()}.vep.tsv
		"""

	else
		"""
        editor="${vcf.simpleName.tokenize('_')[3]}"
        awk '{if (\$4 ~ /^custom/) {print > "custom_regions.bed"} else {print \$4 > "protein.list"}}' ${bed}
        export PERL5LIB=/opt/vep/.vep/Plugins/:\${PERL5LIB:-}
	bcftools sort ${vcf} -Oz -o ${vcf.getSimpleName()}.sorted.vcf.gz
        tabix -p vcf ${vcf.getSimpleName()}.sorted.vcf.gz

        if [ -e protein.list ] ; then
        
                vep -o stdout -i ${vcf.getSimpleName()}.sorted.vcf.gz --flag_pick_allele_gene --cache --offline --species ${species} ${cache_version} --format vcf --tab --force_overwrite --dir_cache /opt/vep/.vep/ --fasta /opt/vep/.genome/${params.genome}.fa --gff /opt/vep/.genome/${params.genome}.gff3.gz --variant_class --nearest gene --gene_phenotype --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory --canonical --protein --biotype  --pubmed --shift_hgvs 0 --allele_number --buffer_size 10000 --custom ${vcf.getSimpleName()}.sorted.vcf.gz,\$editor,vcf,exact,0,Protospacer,PAM,Nchange | filter_vep --filter "SYMBOL in protein.list" -o "${vcf.getSimpleName()}_Protein.vep"
else
                touch "${vcf.getSimpleName()}_Protein.vep"
        fi
        if [ -e custom_regions.bed ] ; then
                bcftools view ${vcf.getSimpleName()}.sorted.vcf.gz -R custom_regions.bed | vep -o "${vcf.getSimpleName()}_region.vep" ${cache_version} --flag_pick --cache --offline --species ${species} --format vcf --tab --force_overwrite --dir_cache /opt/vep/.vep/ --fasta /opt/vep/.genome/${params.genome}.fa  --variant_class --nearest gene --gene_phenotype --ccds --uniprot --hgvs --symbol --numbers --domains --regulatory        
else
                touch ${vcf.getSimpleName()}_region.vep
        fi
        cat "${vcf.getSimpleName()}_Protein.vep" "${vcf.getSimpleName()}_region.vep"  | awk '/^##/ {next} /^#/ && c++ {next} 1' > tmp.tsv
        awk -v value=\$editor 'BEGIN {OFS="\t"} NR==1 {\$0 = \$0 OFS "editor"} NR>1 {\$0 = \$0 OFS value} 1' tmp.tsv > ${vcf.getSimpleName()}.vep.tsv

		"""

}

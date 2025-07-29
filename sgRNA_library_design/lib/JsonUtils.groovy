// lib/JsonUtils.groovy

import groovy.json.JsonSlurper

class JsonUtils {

    static Map loadAssemblyParams(String jsonFile, String assembly) {
        def json        = new JsonSlurper().parse(new File(jsonFile))
        if (!assembly) assembly = json.keySet().first()

        def entry = json[assembly] 
                ?: json.find { k,v -> v.aliases?.contains(assembly) }?.value

        if (!entry) {
            def available = json.collect { k,v ->
                ([k] + (v.aliases ?: [])).join(', ')
            }.join(' | ')
            throw new IllegalArgumentException(
                    "Invalid genome '${assembly}'. Available: ${available}")
        }

        return [
            genome            : assembly,
            release           : entry.release ?: '',
            source            : entry.source ?: '',
            mane_select       : entry.mane_select ?: '',
            name              : entry.species_info?.name ?: '',
            name_url              : entry.species_info?.url_name ?: '',
            vep_cache_url     : entry.urls?.vep_cache ?: '',
            fasta_url         : entry.urls?.fasta ?: '',
            gff3_url          : entry.urls?.gff3 ?: '',
            Boyle_Lab_url     : entry.urls?.blacklist_regions_url?.Boyle_Lab_url ?: '',
            Dust_repeat_url   : entry.urls?.blacklist_regions_url?.Dust_repeat_url ?: '',
            repeatmasquer_url : entry.urls?.blacklist_regions_url?.repeatmasquer_url ?: '',
            fasta_loc         : entry.paths?.fasta_loc ?: '',
            gffdb_loc         : entry.paths?.gffdb_loc ?: '',
            gff3_loc          : entry.paths?.gff3_loc ?: '',
            vep_cache_loc     : entry.paths?.vep_cache_loc ?: '',
            bowtie_index_loc  : entry.paths?.bowtie_index_loc ?: '',
            bsgenome_loc      : entry.paths?.bsgenome_loc ?: '',
            Boyle_Lab         : entry.paths?.blacklist_regions?.Boyle_Lab ?: '',
            Dust_repeat       : entry.paths?.blacklist_regions?.Dust_repeat ?: '',
            repeatmasker      : entry.paths?.blacklist_regions?.repeatmasquer ?: ''
        ]
    }

    static Map loadCasVariant(String casVariantJsonPath, String casName) {
        def casDb = new JsonSlurper().parse(new File(casVariantJsonPath))
        def entry = casDb[casName]
        if (!entry) {
            def available = casDb.keySet().join(', ')
            throw new IllegalArgumentException(
                    "Invalid casName '${casName}'. Available: ${available}")
        }
        return [
            pam        : entry.pam,
            length     : entry.length,
            pamLength  : entry.pamLength,
            CFD_access : entry.CFD_access
        ]
    }
}


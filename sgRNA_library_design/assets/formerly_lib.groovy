    static boolean checkGenomeJsonAssembly(String jsonFile, String assembly) {
        def jsonSlurper = new JsonSlurper()
        def file       = new File(jsonFile)
        def exists     = true

        if (!file.exists()) {
            println "✗ File does not exist: ${jsonFile}"
            exists = false
        } else {
            def json = jsonSlurper.parse(file)
            println "Checking paths for ${assembly} using ${jsonFile}"
            println "\nTop-level JSON keys:"
            json.keySet().each { println it }

            println "\nKey-existence check for '${assembly}':"
            if (!json.containsKey(assembly)) {
                println "✗ MISSING key '${assembly}'"
                exists = false
            } else {
                println "✓ Found key '${assembly}'"
                def entry        = json[assembly]
                def requiredKeys = ['fasta_loc','gffdb_loc','bsgenome_loc','bowtie_index_loc']
                def optionalKeys = ['Boyle_Lab','Dust_repeat','repeatmasquer','vep_cache_loc']

                requiredKeys.each { key ->
                    def path = entry.paths[key]
                    if (!path) {
                        println "Missing required key or empty value: ${key}"
                        exists = false
                    } else if (!new File(path).exists()) {
                        println "Required ${key} file not found: ${path}"
                        exists = false
                    } else {
                        println "${key} file found at: ${path}"
                    }
                }

                optionalKeys.each { key ->
                    def path = key=='vep_cache_loc'
                               ? entry.paths[key]
                               : entry.paths.blacklist_regions[key]
                    if (!path) {
                        println "[Note] Optional '${key}' not provided."
                    } else if (!new File(path).exists()) {
                        println "[Note] '${key}' set to '${path}' but file does not exist."
                        exists = false
                    } else {
                        println "[Note] Optional '${key}' file found at: ${path}"
                    }
                }
            }
        }

        return exists
    }

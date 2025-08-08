// lib/ParamUtils.groovy

class ParamUtils {
static Closure generateParamHelpMessage = { required, optional, descs ->
    def msg = new StringBuilder("PARAMETER HELP:\n")
        if (required) {
                msg << "\nRequired:\n"
                required.each {
                        msg << "  --${it}" + (descs[it] ? " : ${descs[it]}" : "") + "\n"
                }
        }
        if (optional) {
                msg << "\nOptional:\n"
                optional.each {
                        msg << "  --${it}" + (descs[it] ? " : ${descs[it]}" : "") + "\n"
                }
        }

        return msg.toString()
        }
static void validateExactParamSet(
    Map params,
    List requiredKeys = [],
    List optionalKeys = [],
    Closure helpMessage = null,
    Map paramDescriptions = [:],
    String pipelineDescription = null,
    List ignoredKeys = []
) {
    def allowedKeys = requiredKeys + optionalKeys + ignoredKeys
    def provided = params.keySet().findAll { params[it] != null }
    def illegal = provided - allowedKeys

    def prefix = pipelineDescription ? "↪ ${pipelineDescription}\n\n" : ""

    if (illegal) {
        def msg = "${prefix}Unexpected parameter(s): ${illegal.join(', ')}. Allowed: ${requiredKeys + optionalKeys}"
        msg += helpMessage ? "\n\n" + helpMessage(requiredKeys, optionalKeys, paramDescriptions) : ""
        throw new IllegalArgumentException(msg)
    }

    def missing = requiredKeys.findAll { !params.containsKey(it) || params[it] == null }
    if (missing) {
        def msg = "${prefix}Missing required parameter(s): ${missing.join(', ')}"
        msg += helpMessage ? "\n\n" + helpMessage(requiredKeys, optionalKeys, paramDescriptions) : ""
        throw new IllegalArgumentException(msg)
    }
}
}

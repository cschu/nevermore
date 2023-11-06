

def get_single_input_dir {
	if (params.input.local_dir && params.input.remote_dir) {
		log.info """
			Cannot process both --input.local_dir and --input.remote_dir. Please check input parameters.
		""".stripIndent()
		exit 1
	} else if (!params.input.local_dir && !params.input.remote_dir) {
		log.info """
			Please set either --input.local_dir or --input.remote_dir.
		""".stripIndent()
		exit 1
	}
	return (params.input_dir) ? params.input_dir : params.remote_input_dir
}

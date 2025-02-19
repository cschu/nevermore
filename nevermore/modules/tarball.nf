// EXPERIMENTAL

process generate_tarball_old {

	input:
	path(files)

	output:
	path("results.tar.gz")

	script:
	"""
	mkdir results
    ls ${files} | xargs -I {} sh -c 'd=\$(echo {} | sed "s/\\..\\+\$//") && mkdir -p results/\$d && mv -v {} results/\$d/'

    tar chvzf results.tar.gz results
	"""


}

params.restype = "profiles"

process generate_tarball {

    input:
    path("results/${params.restype}/")
    output:
    path("results.tar.gz")

    script:
    """
    tar chvzf results.tar.gz results
    """

}
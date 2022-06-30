import os
import sys
import argparse

walker = os.walk(sys.argv[1])
print("sample", "r1", "r2", "singles", sep=",")
for cwd, dirs, files in walker:
	fastqs = [f for f in files if f.endswith(".fastq.gz")]
	if fastqs:
		sample = os.path.basename(cwd)
		fastqs = sorted(fastqs)
		if 'single' in fastqs[0] and len(fastqs) > 1:
			fastqs.append(fastqs.pop(0))
		fastqs = [os.path.abspath(os.path.join(cwd, f)) for f in fastqs] + ["", ""]
		print(sample, *fastqs[:3], sep=",")

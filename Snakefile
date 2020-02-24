import glob
import pandas as pandas

from snakemake.utils import min_version

#### GLOBAL PARAMETERS AND CONFIG ####
min_version('5.4')
configfile: "config.yaml"



#### GLOBAL SCOPE FUNCTIONS ####
def get_resource(rule,resource):
	'''
	Attempt to parse config.yaml to retrieve resources available for each rule.
	It will revert to default if a key error is found
	'''
	try:
		return config['rules'][rule]['res'][resource]
	except KeyError:
		print(f'Failed to resolve resource config for rule {rule}/{resource}: default parameters have been set')
		return config['rules']['default']['res'][resource]




#	
#   
#
rule filter_ccle_mutations:
	input:
		config['ccle_datasets'] + '/CCLE_mutations.csv'
	output:
		"filtered_mutations.csv"
	script:
	
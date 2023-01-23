### The first script here adjust the microarray annotation file.
#################################### Python part to adjust annotation
####################################
# when a probe points to multiple genes I concatenate the gene names into a single string.
full_file = open('MoGene2_081020W_MB01-06_with annotation.RMA-GENE-FULL - Group 1.TXT', 'r')
translated_file = open('MoGene2_081020W_MB01-06_with_proper_annotation.tsv', 'w')
#translated_file = open('MoGene2_081020W_annotation.tsv', 'w')
for line in full_file:
#	splitted = line[:-1].split('\t')[0:9]
	splitted = line[:-1].split('\t')[1:9]
	gene_names = splitted[7].split('//')
#	gene_names = splitted[8].split('//')
	gene_names_set = set()
	for i in gene_names:
		gene_names_set.add(i.strip())  #remove duplicates  ##.strip() removes leading and trailing whitespaces
	gene_name = '/'.join([a for a in gene_names_set])  #properly describe probes that represent more than 1 gene
	if(gene_name != '---'):
#		dontprint = translated_file.write(splitted[0] + '\t' +  gene_name + '\n')
		dontprint = translated_file.write(gene_name + '\t' + '\t'.join(splitted[0:6]) + '\n')

full_file.close()
translated_file.close()

print('done, now run script 002')

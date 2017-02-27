data=open('1','r')
result=open('1_modified.chain','w')

roman2arabic = {'chrI':'chr1','chrII':'chr2','chrIII':'chr3','chrIV':'chr4','chrV':'chr5','chrVI':'chr6'
,'chrVII':'chr7','chrVIII':'chr8','chrIX':'chr9','chrX':'chr10','chrXI':'chr11','chrXII':'chr12','chrXIII':'chr13'
,'chrXIV':'chr14','chrXV':'chr15','chrXVI':'chr16',}


for records in data:
	tmp=records.split(' ')
	if len(tmp)!=0:
		if tmp[0]=='chain':
			for rec in tmp:
				if rec in roman2arabic:
					result.write(roman2arabic[rec]+' ')
				else:
					result.write(rec+' ')
		else:
			result.write(records)
	else:
		result.write(records)

result.close()
data.close()









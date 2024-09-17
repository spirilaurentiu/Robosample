from modeller import *
from modeller.automodel import *

# create ali file
seq = "GPLGSPEFAQNITARIGEPLVLKCKGAPKKPPQRLEWKLNTGRTEAWKVLSPQGGGPWDSVARVLPNGSLFLPAVGIQDEGIFRCQAMNRNGKETKSNYRVRVYQIPGKPEIVDSASELTAGVPNKVGTCVSEGSYPAGTLSWHLDGKPLVPNEKGVSVKEQTRRHPETGLFTLQSELMVTPARGGDPRPTFSCSFSPGLPRHRALRTAPIQPRVWEPHHHHHH"
seq += '*'
name = 'rage'
template = '3cjj'

# write ali file
ali_file = f'{name}.ali'
with open(ali_file, 'w') as f:
    f.write(f">P1;{name}\n")
    f.write(f"sequence:{name}:::::::0.00: 0.00\n")
    for i in range(0, len(seq), 75):
        f.write(seq[i:i+75] + '\n')

env = Environ()
aln = Alignment(env)
mdl = Model(env, file=template, model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes=f'{template}A', atom_files=f'{template}.pdb')
aln.append(file=ali_file, align_codes=name)
aln.align2d(max_gap_length=50)

# align sequence to template
aln.write(file=f'{name}-{template}.ali', alignment_format='PIR')

# build model
a = AllHModel(env, alnfile=f'{name}-{template}.ali',
              knowns=f'{template}A', sequence=name,
              assess_methods=(assess.DOPE,
                              assess.GA341))
a.starting_model = 1
a.ending_model = 1
a.make()

# get name of the model
model_name = f'{name}.B99990001.pdb'
print(a.name)

# protonate

# ace nme cap
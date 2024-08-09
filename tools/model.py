from modeller import *
from modeller.automodel import *

env = Environ()
aln = Alignment(env)
mdl = Model(env, file='3cjj', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='3cjjA', atom_files='3cjj.pdb')
aln.append(file='rage.ali', align_codes='rage')
aln.align2d(max_gap_length=50)
aln.write(file='rage-3cjj.ali', alignment_format='PIR')
aln.write(file='rage-3cjj.pap', alignment_format='PAP')

a = AllHModel(env, alnfile='rage-3cjj.ali',
              knowns='3cjjA', sequence='rage',
              assess_methods=(assess.DOPE,
                              assess.GA341))
a.starting_model = 1
a.ending_model = 1
a.make()

# protonate

# ace nme cap